package motif

import (
	"context"
	"errors"
	"fmt"
	"os"
	"sync"
)

// TaskStatus 描述任务状态。
type TaskStatus string

const (
	StatusRunning   TaskStatus = "running"
	StatusCompleted TaskStatus = "completed"
	StatusCanceled  TaskStatus = "canceled"
	StatusFailed    TaskStatus = "failed"
)

// ErrCanceled 表示任务因显式取消而在工作单元边界中断。
var ErrCanceled = errors.New("motif: task canceled")

// Snapshot 是任务状态的只读快照。
type Snapshot struct {
	ID       string
	Status   TaskStatus
	Progress float64
	Result   *Result
	Error    string
}

type taskState struct {
	mu       sync.Mutex
	status   TaskStatus
	progress float64
	err      string
	result   *Result

	// cancel 在工作单元边界被关闭以请求中断。
	cancel chan struct{}

	// done 在任务 goroutine 结束时关闭，供 Wait 阻塞。
	done chan struct{}
}

// Service 管理多个可并发、可中断恢复的发现任务。
// 任务状态持久化在 stateDir 中：每个任务一个 JSON 检查点文件，
// 原子写入保证中断时不会出现部分结果。
type Service struct {
	dir   string
	mu    sync.Mutex
	tasks map[string]*taskState
}

// NewService 创建基于 stateDir 持久化检查点的服务。
func NewService(stateDir string) (*Service, error) {
	if stateDir == "" {
		return nil, fmt.Errorf("service: stateDir must not be empty")
	}
	if err := os.MkdirAll(stateDir, 0o755); err != nil {
		return nil, fmt.Errorf("service: cannot create state dir: %w", err)
	}
	return &Service{dir: stateDir, tasks: make(map[string]*taskState)}, nil
}

// Start 启动一个新的发现任务。
// 相同 ID 的任务（无论已完成还是已中断）不能重复创建，应使用 Resume。
func (s *Service) Start(ctx context.Context, id string, seqs []Sequence, cfg Config) error {
	if id == "" {
		return fmt.Errorf("service: task id must not be empty")
	}
	s.mu.Lock()
	if _, exists := s.tasks[id]; exists {
		s.mu.Unlock()
		return fmt.Errorf("service: task %q already exists; use Resume to continue it", id)
	}
	if _, err := loadCheckpoint(s.dir, id); err == nil {
		s.mu.Unlock()
		return fmt.Errorf("service: task %q has persisted state; use Resume to continue it", id)
	} else if !errors.Is(err, os.ErrNotExist) {
		s.mu.Unlock()
		return fmt.Errorf("service: cannot inspect task %q: %w", id, err)
	}

	state, err := newDiscoveryState(seqs, cfg)
	if err != nil {
		s.mu.Unlock()
		return err
	}
	cp := &checkpoint{
		ID:         id,
		Status:     string(StatusRunning),
		Config:     state.cfg,
		Sequences:  seqs,
		NullScores: [][]float64{},
	}
	if err := cp.store(s.dir); err != nil {
		s.mu.Unlock()
		return fmt.Errorf("service: cannot write checkpoint: %w", err)
	}

	ts := &taskState{status: StatusRunning, cancel: make(chan struct{}), done: make(chan struct{})}
	s.tasks[id] = ts
	s.mu.Unlock()

	go s.runLoop(ctx, id, ts, state)
	return nil
}

// Resume 从中断点恢复任务并继续运行直到完成。
// 恢复的任务重放检查点：已完成的工作单元直接复用，未完成单元逐位重演，
// 最终输出与不中断运行完全一致。
func (s *Service) Resume(ctx context.Context, id string) error {
	s.mu.Lock()
	var existing *taskState
	if ts, exists := s.tasks[id]; exists {
		ts.mu.Lock()
		running := ts.status == StatusRunning
		ts.mu.Unlock()
		if running {
			s.mu.Unlock()
			return fmt.Errorf("service: task %q is already running", id)
		}
		if ts.status == StatusCompleted {
			s.mu.Unlock()
			return fmt.Errorf("service: task %q already completed", id)
		}
		// 进程内已取消/失败的任务：用新 goroutine 替换旧条目后继续。
		existing = ts
	}

	cp, err := loadCheckpoint(s.dir, id)
	if err != nil {
		s.mu.Unlock()
		if errors.Is(err, os.ErrNotExist) {
			return fmt.Errorf("service: unknown task %q", id)
		}
		return err
	}
	if TaskStatus(cp.Status) == StatusCompleted {
		s.mu.Unlock()
		return fmt.Errorf("service: task %q already completed", id)
	}

	state, err := restoreState(cp)
	if err != nil {
		s.mu.Unlock()
		return err
	}
	ts := &taskState{status: StatusRunning, cancel: make(chan struct{}), done: make(chan struct{})}
	s.tasks[id] = ts
	s.mu.Unlock()
	if existing != nil {
		// 等待旧 goroutine 完全退出后再启动新的，避免并发写检查点。
		<-existing.done
	}

	go s.runLoop(ctx, id, ts, state)
	return nil
}

// Cancel 请求中断任务。中断只在工作单元边界生效：
// 当前单元完成并落盘后 goroutine 退出，不会留下部分结果。
func (s *Service) Cancel(id string) error {
	s.mu.Lock()
	ts, ok := s.tasks[id]
	s.mu.Unlock()
	if !ok {
		if cp, err := loadCheckpoint(s.dir, id); err == nil && TaskStatus(cp.Status) != StatusCompleted {
			cp.Status = string(StatusCanceled)
			return cp.store(s.dir)
		}
		return fmt.Errorf("service: unknown running task %q", id)
	}
	ts.mu.Lock()
	if ts.status != StatusRunning {
		status := ts.status
		ts.mu.Unlock()
		return fmt.Errorf("service: task %q not running (status %q)", id, status)
	}
	close(ts.cancel)
	ts.mu.Unlock()
	return nil
}

// Get 返回任务当前状态快照。
func (s *Service) Get(id string) (*Snapshot, error) {
	s.mu.Lock()
	ts, ok := s.tasks[id]
	s.mu.Unlock()
	if ok {
		ts.mu.Lock()
		defer ts.mu.Unlock()
		return &Snapshot{
			ID:       id,
			Status:   ts.status,
			Progress: ts.progress,
			Result:   ts.result,
			Error:    ts.err,
		}, nil
	}

	cp, err := loadCheckpoint(s.dir, id)
	if err != nil {
		return nil, fmt.Errorf("service: unknown task %q", id)
	}
	snap := &Snapshot{ID: id, Status: TaskStatus(cp.Status), Error: cp.Error}
	if cp.Result != nil {
		snap.Result = fromSavedResult(*cp.Result)
		snap.Progress = 1
	} else {
		snap.Progress = checkpointProgress(cp)
	}
	return snap, nil
}

// Wait 阻塞直到任务完成、失败或 ctx 取消。
func (s *Service) Wait(ctx context.Context, id string) (*Snapshot, error) {
	s.mu.Lock()
	ts, ok := s.tasks[id]
	s.mu.Unlock()
	if !ok {
		return s.Get(id)
	}
	select {
	case <-ts.done:
		return s.Get(id)
	case <-ctx.Done():
		return nil, ctx.Err()
	}
}

// runLoop 驱动一个任务状态机，在每个单元边界落盘并响应取消。
func (s *Service) runLoop(ctx context.Context, id string, ts *taskState, state *discoveryState) {
	defer close(ts.done)

	hook := func(phase string, done, total int) error {
		progress := unitProgress(phase, done, total)
		ts.mu.Lock()
		ts.progress = progress
		ts.mu.Unlock()

		cp := s.buildCheckpoint(id, state)
		if err := cp.store(s.dir); err != nil {
			return fmt.Errorf("persist checkpoint: %w", err)
		}

		select {
		case <-ts.cancel:
			return ErrCanceled
		case <-ctx.Done():
			return ctx.Err()
		default:
			return nil
		}
	}

	result, err := state.run(hook)

	ts.mu.Lock()
	defer ts.mu.Unlock()
	switch {
	case errors.Is(err, ErrCanceled):
		ts.status = StatusCanceled
		cp := s.buildCheckpoint(id, state)
		cp.Status = string(StatusCanceled)
		_ = cp.store(s.dir)
	case err != nil:
		ts.status = StatusFailed
		ts.err = err.Error()
		cp := s.buildCheckpoint(id, state)
		cp.Status = string(StatusFailed)
		cp.Error = err.Error()
		_ = cp.store(s.dir)
	default:
		ts.status = StatusCompleted
		ts.progress = 1
		ts.result = result
		cp := s.buildCheckpoint(id, state)
		cp.Status = string(StatusCompleted)
		sr := toSavedResult(result)
		cp.Result = &sr
		if err := cp.store(s.dir); err != nil {
			ts.status = StatusFailed
			ts.err = err.Error()
		}
	}
}

func (s *Service) buildCheckpoint(id string, state *discoveryState) *checkpoint {
	saved := make([]savedModel, 0, len(state.models))
	for _, m := range state.models {
		if m != nil {
			saved = append(saved, toSavedModel(m))
		}
	}
	nullScores := make([][]float64, len(state.models))
	for k := range state.models {
		if k < len(state.nullScores) {
			nullScores[k] = append([]float64(nil), state.nullScores[k]...)
		} else {
			nullScores[k] = []float64{}
		}
	}
	return &checkpoint{
		ID:               id,
		Status:           string(StatusRunning),
		Config:           state.cfg,
		Sequences:        state.seqs,
		DiscoverDone:     state.discoverDone,
		DiscoverComplete: state.discoverComplete,
		SavedModels:      saved,
		ShuffleDone:      state.shuffleDone,
		NullScores:       nullScores,
	}
}

// restoreState 从检查点重建内存中的发现状态。
func restoreState(cp *checkpoint) (*discoveryState, error) {
	state, err := newDiscoveryState(cp.Sequences, cp.Config)
	if err != nil {
		return nil, err
	}
	state.discoverDone = cp.DiscoverDone
	state.discoverComplete = cp.DiscoverComplete
	state.shuffleDone = cp.ShuffleDone
	if state.discoverComplete {
		models := make([]*emModel, 0, len(cp.SavedModels))
		for _, sm := range cp.SavedModels {
			models = append(models, fromSavedModel(sm))
		}
		state.models = models
	} else {
		models := make([]*emModel, cp.Config.ExpectMotifs)
		for i, sm := range cp.SavedModels {
			if i < len(models) {
				models[i] = fromSavedModel(sm)
			}
		}
		state.models = models
	}
	if cp.NullScores != nil {
		state.nullScores = cp.NullScores
	}
	if cp.Result != nil {
		state.result = fromSavedResult(*cp.Result)
	}
	return state, nil
}

// unitProgress 将两阶段进度映射到 [0,1)：发现阶段占 30%，重排阶段占 70%。
func unitProgress(phase string, done, total int) float64 {
	switch phase {
	case "discover":
		return 0.3 * float64(done) / float64(total)
	case "shuffle":
		return 0.3 + 0.7*float64(done)/float64(total)
	default:
		return 0
	}
}

func checkpointProgress(cp *checkpoint) float64 {
	d := 0.3 * float64(cp.DiscoverDone) / float64(cp.Config.ExpectMotifs)
	if cp.Config.Shuffles > 0 {
		d += 0.7 * float64(cp.ShuffleDone) / float64(cp.Config.Shuffles)
	}
	if d > 1 {
		d = 1
	}
	return d
}
