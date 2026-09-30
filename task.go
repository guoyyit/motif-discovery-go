package motif

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sync"
	"time"
)

// TaskStatus 表示发现任务的生命周期状态。
type TaskStatus int

const (
	StatusPending TaskStatus = iota
	StatusRunning
	StatusPaused
	StatusDone
	StatusFailed
	StatusCanceled
)

// String 返回状态的可读名称。
func (s TaskStatus) String() string {
	switch s {
	case StatusPending:
		return "pending"
	case StatusRunning:
		return "running"
	case StatusPaused:
		return "paused"
	case StatusDone:
		return "done"
	case StatusFailed:
		return "failed"
	case StatusCanceled:
		return "canceled"
	default:
		return "unknown"
	}
}

// ErrPaused 是任务在检查点边界被中断时 run 循环返回的哨兵错误。
var ErrPaused = errors.New("motif: 任务已在检查点暂停")

// ErrCanceled 是任务被取消时 run 循环返回的哨兵错误。
var ErrCanceled = errors.New("motif: 任务已取消")

type roundRecord struct {
	Starts []int   `json:"starts"`
	Score  float64 `json:"score"`
}

// Checkpoint 是发现任务可持久化的断点状态。
// 每完成一个初始化轮次或一次随机重排即更新检查点，保证：
//   - 中断/恢复不产生部分结果（结果只在全部完成时一次性落盘）；
//   - 恢复不重复计算（已完成的轮次/重放编号被记录并跳过）；
//   - 相同种子与输入在是否中断两种情况下输出逐位一致。
type Checkpoint struct {
	mu sync.Mutex `json:"-"`

	TaskID     string     `json:"task_id"`
	Seed       uint64     `json:"seed"`
	ConfigJSON string     `json:"config_json"`
	Status     TaskStatus `json:"status"`
	UpdatedAt  time.Time  `json:"updated_at"`

	DoneRounds  map[string]bool        `json:"done_rounds"`
	SlotBest    map[string]roundRecord `json:"slot_best"`
	Motifs      []Motif                `json:"motifs"`
	Blocked     []map[int]bool         `json:"blocked"`
	NullSamples map[string][]float64   `json:"null_samples"`

	pauseRequested bool    `json:"-"`
	canceled       bool    `json:"-"`
	result         *Result `json:"-"`
	failure        error   `json:"-"`
	path           string  `json:"-"`
}

func newCheckpoint(cfg Config) *Checkpoint {
	digest := configDigest(cfg)
	return &Checkpoint{
		TaskID:      digest[:16],
		Seed:        cfg.Seed,
		ConfigJSON:  digest,
		Status:      StatusPending,
		DoneRounds:  map[string]bool{},
		SlotBest:    map[string]roundRecord{},
		NullSamples: map[string][]float64{},
	}
}

func configDigest(cfg Config) string {
	b, _ := json.Marshal(cfg)
	sum := sha256.Sum256(b)
	return hex.EncodeToString(sum[:])
}

func (c *Checkpoint) match(cfg Config) error {
	if c.Seed != cfg.Seed {
		return fmt.Errorf("motif: 检查点种子 %d 与配置种子 %d 不一致", c.Seed, cfg.Seed)
	}
	if c.ConfigJSON != configDigest(cfg) {
		return errors.New("motif: 检查点配置与当前配置不一致，拒绝恢复")
	}
	return nil
}

func slotRoundKey(slot, round int) string {
	return fmt.Sprintf("%d:%d", slot, round)
}

func (c *Checkpoint) isRoundDone(slot, round int) bool {
	c.mu.Lock()
	defer c.mu.Unlock()
	return c.DoneRounds[slotRoundKey(slot, round)]
}

func (c *Checkpoint) recordRound(slot, round, nSeq int, m *Motif) {
	c.mu.Lock()
	defer c.mu.Unlock()
	key := slotRoundKey(slot, round)
	c.DoneRounds[key] = true
	starts := make([]int, nSeq)
	for _, occ := range m.Sites {
		starts[occ.SeqIndex] = occ.Start
	}
	rec := roundRecord{Starts: starts, Score: m.Score}
	bestKey := fmt.Sprintf("best:%d", slot)
	prev, ok := c.SlotBest[bestKey]
	if !ok || m.Score > prev.Score {
		c.SlotBest[bestKey] = rec
	}
	c.touchLocked()
	_ = c.persistLocked()
}

// bestStarts 返回某槽位当前最优轮次的位点分配，供恢复时重建 PWM。
func (c *Checkpoint) bestStarts(slot int) ([]int, bool) {
	c.mu.Lock()
	defer c.mu.Unlock()
	rec, ok := c.SlotBest[fmt.Sprintf("best:%d", slot)]
	if !ok {
		return nil, false
	}
	return append([]int(nil), rec.Starts...), true
}

func (c *Checkpoint) recordMotif(m Motif, blocked []map[int]bool) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.Motifs = append(c.Motifs, m)
	c.Blocked = cloneBlocked(blocked)
	c.touchLocked()
	_ = c.persistLocked()
}

func (c *Checkpoint) recordMotifUpdate(slot int, m Motif) {
	c.mu.Lock()
	defer c.mu.Unlock()
	if slot < len(c.Motifs) {
		c.Motifs[slot] = m
	}
	c.touchLocked()
	_ = c.persistLocked()
}

func (c *Checkpoint) nullStats(slot int) []float64 {
	c.mu.Lock()
	defer c.mu.Unlock()
	return append([]float64(nil), c.NullSamples[fmt.Sprintf("null:%d", slot)]...)
}

func (c *Checkpoint) recordNull(slot int, stats []float64) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.NullSamples[fmt.Sprintf("null:%d", slot)] = append([]float64(nil), stats...)
	c.touchLocked()
	_ = c.persistLocked()
}

func (c *Checkpoint) checkCancel() error {
	c.mu.Lock()
	defer c.mu.Unlock()
	if c.canceled {
		return ErrCanceled
	}
	if c.pauseRequested {
		c.Status = StatusPaused
		c.touchLocked()
		_ = c.persistLocked()
		return ErrPaused
	}
	return nil
}

func (c *Checkpoint) complete(res *Result) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.result = res
	c.Status = StatusDone
	c.touchLocked()
	_ = c.persistLocked()
}

func (c *Checkpoint) touchLocked() {
	c.UpdatedAt = time.Now().UTC()
}

// persistLocked 将检查点原子写入 path（先写临时文件再 rename）。
func (c *Checkpoint) persistLocked() error {
	if c.path == "" {
		return nil
	}
	b, err := json.Marshal(c)
	if err != nil {
		return err
	}
	tmp := c.path + ".tmp"
	if err := os.WriteFile(tmp, b, 0o600); err != nil {
		return err
	}
	return os.Rename(tmp, c.path)
}

// Save 显式将检查点持久化到 dir/<TaskID>.json（原子写）。
func (c *Checkpoint) Save(dir string) (string, error) {
	c.mu.Lock()
	defer c.mu.Unlock()
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return "", err
	}
	c.path = filepath.Join(dir, c.TaskID+".json")
	if err := c.persistLocked(); err != nil {
		return "", err
	}
	return c.path, nil
}

// LoadCheckpoint 从 JSON 文件读取检查点。
func LoadCheckpoint(path string) (*Checkpoint, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	c := newCheckpoint(DefaultConfig())
	if err := json.Unmarshal(b, c); err != nil {
		return nil, err
	}
	if c.DoneRounds == nil {
		c.DoneRounds = map[string]bool{}
	}
	if c.SlotBest == nil {
		c.SlotBest = map[string]roundRecord{}
	}
	if c.NullSamples == nil {
		c.NullSamples = map[string][]float64{}
	}
	c.path = path
	c.pauseRequested = false
	return c, nil
}

func cloneBlocked(in []map[int]bool) []map[int]bool {
	out := make([]map[int]bool, len(in))
	for i := range in {
		out[i] = make(map[int]bool, len(in[i]))
		for k, v := range in[i] {
			out[i][k] = v
		}
	}
	return out
}
