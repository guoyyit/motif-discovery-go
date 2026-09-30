package motif

import (
	"context"
	"sync"
)

// Task 表示一个可中断、可恢复、可取消的发现任务。
type Task struct {
	cp   *Checkpoint
	seqs []Seq
	cfg  Config
	mu   sync.Mutex
	wg   sync.WaitGroup
}

// TaskManager 管理多个并发运行、彼此隔离的发现任务。
// 每个任务持有独立的输入、配置与检查点；随机性完全来自确定性 RNG，
// 因此同一任务重跑逐位一致，不同任务之间互不干扰。
type TaskManager struct {
	mu    sync.Mutex
	tasks map[string]*Task
}

// NewTaskManager 构造空的任务管理器。
func NewTaskManager() *TaskManager {
	return &TaskManager{tasks: map[string]*Task{}}
}

// SubmitWithCheckpoint 提交使用外部检查点（含磁盘路径）的任务，
// 检查点在每个边界原子落盘，可用于跨进程中断/恢复。
func (m *TaskManager) SubmitWithCheckpoint(seqs []Seq, cfg Config, cp *Checkpoint) (string, error) {
	if err := cp.match(cfg); err != nil {
		return "", err
	}
	task := &Task{cp: cp, seqs: cloneSeqs(seqs), cfg: cfg}
	m.mu.Lock()
	m.tasks[cp.TaskID] = task
	m.mu.Unlock()
	task.wg.Add(1)
	go task.run()
	return cp.TaskID, nil
}

// Submit 提交一个新的发现任务并异步执行，返回任务 ID（即配置+种子的摘要前缀）。
func (m *TaskManager) Submit(seqs []Seq, cfg Config) (string, error) {
	if err := cfg.Validate(); err != nil {
		return "", err
	}
	cp := newCheckpoint(cfg)
	seqsCopy := cloneSeqs(seqs)
	task := &Task{
		cp:   cp,
		seqs: seqsCopy,
		cfg:  cfg,
	}
	m.mu.Lock()
	if _, exists := m.tasks[cp.TaskID]; exists {
		m.mu.Unlock()
		return "", errTaskExists(cp.TaskID)
	}
	m.tasks[cp.TaskID] = task
	m.mu.Unlock()

	task.wg.Add(1)
	go task.run()
	return cp.TaskID, nil
}

// SubmitResumable 从已有检查点恢复任务（配置必须与检查点匹配）。
func (m *TaskManager) SubmitResumable(seqs []Seq, cfg Config, cp *Checkpoint) (string, error) {
	if err := cp.match(cfg); err != nil {
		return "", err
	}
	task := &Task{
		cp:   cp,
		seqs: cloneSeqs(seqs),
		cfg:  cfg,
	}
	cp.Status = StatusRunning
	m.mu.Lock()
	m.tasks[cp.TaskID] = task
	m.mu.Unlock()
	task.wg.Add(1)
	go task.run()
	return cp.TaskID, nil
}

func (t *Task) run() {
	defer t.wg.Done()
	t.cp.mu.Lock()
	t.cp.Status = StatusRunning
	t.cp.mu.Unlock()

	_, err := DiscoverWithProgress(t.seqs, t.cfg, t.cp)

	t.cp.mu.Lock()
	defer t.cp.mu.Unlock()
	switch err {
	case nil:
		t.cp.Status = StatusDone
	case ErrPaused:
		t.cp.Status = StatusPaused
	case ErrCanceled:
		t.cp.Status = StatusCanceled
	default:
		t.cp.Status = StatusFailed
		t.cp.failure = err
	}
}

// Get 返回指定任务的当前状态与（完成后的）结果。
func (m *TaskManager) Get(id string) (TaskStatus, *Result, error) {
	t, err := m.require(id)
	if err != nil {
		return StatusFailed, nil, err
	}
	t.cp.mu.Lock()
	defer t.cp.mu.Unlock()
	return t.cp.Status, t.cp.result, t.cp.failure
}

// Pause 请求在最近的检查点边界中断任务。
func (m *TaskManager) Pause(id string) error {
	t, err := m.require(id)
	if err != nil {
		return err
	}
	t.cp.mu.Lock()
	defer t.cp.mu.Unlock()
	if t.cp.Status != StatusRunning && t.cp.Status != StatusPending {
		return errBadStatus(id, t.cp.Status)
	}
	t.cp.pauseRequested = true
	return nil
}

// Resume 从检查点恢复任务，输出与不中断运行逐位一致。
func (m *TaskManager) Resume(id string) error {
	t, err := m.require(id)
	if err != nil {
		return err
	}
	t.cp.mu.Lock()
	defer t.cp.mu.Unlock()
	if t.cp.Status != StatusPaused {
		return errBadStatus(id, t.cp.Status)
	}
	t.cp.pauseRequested = false
	t.cp.Status = StatusRunning
	t.wg.Add(1)
	go t.run()
	return nil
}

// Cancel 取消任务（幂等）。
func (m *TaskManager) Cancel(_ context.Context, id string) error {
	t, err := m.require(id)
	if err != nil {
		return err
	}
	t.cp.mu.Lock()
	defer t.cp.mu.Unlock()
	t.cp.canceled = true
	t.cp.pauseRequested = false
	return nil
}

// Wait 阻塞直到任务结束（完成/失败/取消）或 ctx 结束。
func (m *TaskManager) Wait(ctx context.Context, id string) (TaskStatus, *Result, error) {
	t, err := m.require(id)
	if err != nil {
		return StatusFailed, nil, err
	}
	t.cp.mu.Lock()
	status := t.cp.Status
	t.cp.mu.Unlock()
	waitCh := make(chan struct{})
	go func() {
		t.wg.Wait()
		close(waitCh)
	}()
	select {
	case <-waitCh:
	case <-ctx.Done():
		return status, nil, ctx.Err()
	}
	return m.Get(id)
}

// Checkpoint 暴露任务的检查点，便于持久化或检查状态。
func (m *TaskManager) Checkpoint(id string) (*Checkpoint, error) {
	t, err := m.require(id)
	if err != nil {
		return nil, err
	}
	return t.cp, nil
}

func (m *TaskManager) require(id string) (*Task, error) {
	m.mu.Lock()
	defer m.mu.Unlock()
	t, ok := m.tasks[id]
	if !ok {
		return nil, errUnknownTask(id)
	}
	return t, nil
}

func cloneSeqs(in []Seq) []Seq {
	out := make([]Seq, len(in))
	for i, s := range in {
		out[i] = Seq{
			Name: s.Name,
			Desc: s.Desc,
			Data: append([]byte(nil), s.Data...),
		}
	}
	return out
}
