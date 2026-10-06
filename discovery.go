package motif

import (
	"fmt"
	"math"
)

// Config 是模体发现流程的参数配置。
type Config struct {
	// Seed 是所有随机过程的显式主种子，相同种子与输入保证逐位一致的输出。
	Seed uint64
	// MotifLength 是待发现模体的长度 L。
	MotifLength int
	// ExpectMotifs 是期望发现的模体数量 K。
	ExpectMotifs int
	// InitRounds 是每个模体 EM 的随机初始化轮数，取似然最高者。
	InitRounds int
	// ConvergenceTol 是 EM 对数似然的收敛容差（绝对变化量）；缺省（0）填默认值。
	ConvergenceTol float64
	// ScoreThreshold 是位点归一化得分阈值，合法范围 (0,1]。
	// 使用指针以区分“缺省”（nil，填默认值 0.7）与显式给出的值；
	// 显式给出 <=0 或 >1 的值会直接报错，绝不静默替换为默认值。
	// 可用 motif.Float64(0.8) 构造。
	ScoreThreshold *float64
	// SignificanceAlpha 是经多重检验校正后的显著性阈值（FDR），合法范围 (0,1)。
	// 同样以指针区分缺省（nil，填默认值 0.05）与显式值。
	SignificanceAlpha *float64
	// Shuffles 是显著性零分布的重排次数 B，最小经验 p 值为 1/(B+1)。
	Shuffles int
	// MaxEMIterations 是单次 EM 拟合的迭代上限；0 表示使用默认值。
	// 达到上限仍不满足 ConvergenceTol 时，该任务以可区分的
	// ErrEMNotConverged 错误失败，不会挂起，也不会返回半成品结果。
	MaxEMIterations int
	// Pseudocount 是 PWM 的伪计数强度。
	Pseudocount float64
	// Background 为零值（全零）时自动从输入序列估计。
	Background Background
	// MaskLowComplexity 控制低复杂度/简单重复区域掩码：
	//   - nil（缺省）：默认启用保护，AT 串联重复等简单重复不会被报告为显著模体；
	//   - 显式 true：启用掩码，语义与缺省路径完全一致（两条路径结论相同）；
	//   - 显式 false：关闭掩码，由调用方自行承担简单重复被当作候选的风险。
	// 缺省路径与显式开启后的结论逐位一致。
	MaskLowComplexity *bool
}

// Bool 是设置 Config.MaskLowComplexity 等 *bool 字段的辅助函数：
// motif.Bool(true) / motif.Bool(false)。
func Bool(v bool) *bool { return &v }

// Float64 是设置 Config.ScoreThreshold / SignificanceAlpha 等
// *float64 字段的辅助函数：motif.Float64(0.8)。
func Float64(v float64) *float64 { return &v }

// Site 是模体在某条序列中的一个出现位置。
type Site struct {
	SeqIndex int
	Offset   int
	Score    float64
}

// Motif 是单个发现模体及其统计结论。
type Motif struct {
	Index       int
	Matrix      *PWM
	Sites       []Site
	Information float64
	// Statistic 是观测统计量（每条序列最佳归一化得分之和）。
	Statistic float64
	// NullCount 是用于构建经验零分布的重排样本数 B。
	// 结合 PValue 可核对超出次数：#{null >= observed} = PValue*(B+1) - 1。
	NullCount      int
	PValue         float64
	AdjustedPValue float64
	Significant    bool
}

// Result 是一次发现任务的最终输出。
type Result struct {
	Motifs []Motif
}

// 默认参数（零值字段自动填充）。
const (
	defaultMotifLength    = 10
	defaultExpectMotifs   = 3
	defaultInitRounds     = 10
	defaultConvergenceTol = 1e-4
	defaultScoreThreshold = 0.7
	defaultAlpha          = 0.05
	defaultShuffles       = 200
	defaultPseudocount    = 0.5
	// defaultMaxEMIterations 为 EM 单次拟合的默认迭代上限。
	// 缺省配置下拟合始终有界，退化输入也在有限步内结束。
	defaultMaxEMIterations = 1000
)

// maskLowComplexity 返回（填充默认值后的）低复杂度掩码开关的实际取值。
// 缺省（nil）与显式 true 行为一致：均启用保护。
func (c Config) maskLowComplexity() bool {
	return c.MaskLowComplexity == nil || *c.MaskLowComplexity
}

// Validate 检查配置合法性并填充默认值，返回填充后的配置。
//
// 缺省值与非法值严格区分：缺省（整型零值、浮点零值、指针为 nil）视为
// “未配置”，填入默认值；显式给出的越界值（如 ScoreThreshold <= 0、
// SignificanceAlpha 不在 (0,1)、计数参数为负）直接返回错误，错误信息
// 注明字段名与合法取值范围，绝不静默替换为默认值后继续运行。
func (c Config) Validate() (Config, error) {
	// 先填充缺省值（仅限“未配置”语义的零值/nil）。
	if c.MotifLength == 0 {
		c.MotifLength = defaultMotifLength
	}
	if c.ExpectMotifs == 0 {
		c.ExpectMotifs = defaultExpectMotifs
	}
	if c.InitRounds == 0 {
		c.InitRounds = defaultInitRounds
	}
	if c.ConvergenceTol == 0 {
		c.ConvergenceTol = defaultConvergenceTol
	}
	if c.ScoreThreshold == nil {
		c.ScoreThreshold = Float64(defaultScoreThreshold)
	}
	if c.SignificanceAlpha == nil {
		c.SignificanceAlpha = Float64(defaultAlpha)
	}
	if c.Shuffles == 0 {
		c.Shuffles = defaultShuffles
	}
	if c.MaxEMIterations == 0 {
		c.MaxEMIterations = defaultMaxEMIterations
	}
	if c.Pseudocount == 0 {
		c.Pseudocount = defaultPseudocount
	}

	// 再做范围校验：以下值若出现，必然是调用方显式给出的非法值。
	switch {
	case c.MotifLength < 4:
		return c, fmt.Errorf("config: MotifLength must be >= 4 (got %d)", c.MotifLength)
	case c.ExpectMotifs < 1:
		return c, fmt.Errorf("config: ExpectMotifs must be >= 1 (got %d)", c.ExpectMotifs)
	case c.InitRounds < 1:
		return c, fmt.Errorf("config: InitRounds must be >= 1 (got %d)", c.InitRounds)
	case c.Shuffles < 1:
		return c, fmt.Errorf("config: Shuffles must be >= 1 (got %d)", c.Shuffles)
	case c.MaxEMIterations < 1:
		return c, fmt.Errorf("config: MaxEMIterations must be >= 1 (got %d)", c.MaxEMIterations)
	case math.IsNaN(c.ConvergenceTol) || math.IsInf(c.ConvergenceTol, 0):
		return c, fmt.Errorf("config: ConvergenceTol must be a finite number (got %v)", c.ConvergenceTol)
	case c.ConvergenceTol < 0:
		return c, fmt.Errorf("config: ConvergenceTol must be >= 0 (got %v)", c.ConvergenceTol)
	case math.IsNaN(*c.ScoreThreshold) || *c.ScoreThreshold <= 0 || *c.ScoreThreshold > 1:
		return c, fmt.Errorf("config: ScoreThreshold must be within (0,1] (got %v)", *c.ScoreThreshold)
	case math.IsNaN(*c.SignificanceAlpha) || *c.SignificanceAlpha <= 0 || *c.SignificanceAlpha > 1:
		return c, fmt.Errorf("config: SignificanceAlpha must be within (0,1) (got %v)", *c.SignificanceAlpha)
	case math.IsNaN(c.Pseudocount) || math.IsInf(c.Pseudocount, 0):
		return c, fmt.Errorf("config: Pseudocount must be a finite number (got %v)", c.Pseudocount)
	case c.Pseudocount < 0:
		return c, fmt.Errorf("config: Pseudocount must be >= 0 (got %v)", c.Pseudocount)
	}

	// 零值背景保留为零值，由发现流程在运行时从输入序列估计；
	// 显式提供的背景必须是合法概率分布。
	if c.Background != (Background{}) {
		if err := c.Background.Validate(); err != nil {
			return c, err
		}
	}
	return c, nil
}

// progressHook 在每个可恢复工作单元完成时被调用；返回非 nil 表示请求中断。
// phase 描述阶段（"discover"/"shuffle"），done/total 为该阶段内的进度。
type progressHook func(phase string, done, total int) error

// Discover 对序列集合执行模体发现（无中断支持的便捷入口）。
func Discover(seqs []Sequence, cfg Config) (*Result, error) {
	state, err := newDiscoveryState(seqs, cfg)
	if err != nil {
		return nil, err
	}
	return state.run(nil)
}

// discoveryState 是发现流程的可恢复中间状态。
// 每个工作单元（发现模体 k、重排 b）完成后更新该状态；
// 单元之间可安全落盘并中断。每个单元使用独立的确定性随机流，
// 跳过已完成单元不会改变后续单元的随机序列，因此恢复后的输出
// 与不中断运行逐位一致。
type discoveryState struct {
	cfg     Config
	seqs    []Sequence
	working [][]byte

	models           []*emModel
	nullScores       [][]float64
	discoverDone     int
	discoverComplete bool
	shuffleDone      int
	result           *Result
}

func newDiscoveryState(seqs []Sequence, cfg Config) (*discoveryState, error) {
	cfg, err := cfg.Validate()
	if err != nil {
		return nil, err
	}
	if len(seqs) == 0 {
		return nil, fmt.Errorf("discover: no input sequences")
	}

	// 背景基于（掩码前的）输入序列组成估计，重排保持同样组成。
	if cfg.Background == (Background{}) {
		cfg.Background = EstimateBackground(seqs)
	}

	working := make([][]byte, len(seqs))
	for i, seq := range seqs {
		data := append([]byte(nil), seq.Data...)
		if cfg.maskLowComplexity() {
			data = MaskLowComplexity(data)
		}
		working[i] = data
	}

	usable := 0
	for _, data := range working {
		if len(data) >= cfg.MotifLength {
			usable++
		}
	}
	if usable == 0 {
		return nil, fmt.Errorf(
			"discover: no sequence is at least MotifLength (%d) long", cfg.MotifLength)
	}

	return &discoveryState{
		cfg:     cfg,
		seqs:    seqs,
		working: working,
		models:  make([]*emModel, cfg.ExpectMotifs),
	}, nil
}

// run 从当前状态继续执行直到完成；hook 在每个工作单元边界被调用。
func (s *discoveryState) run(hook progressHook) (*Result, error) {
	if !s.discoverComplete {
		if err := s.runDiscover(hook); err != nil {
			return nil, err
		}
	}
	if err := s.runShuffles(hook); err != nil {
		return nil, err
	}
	return s.finalize()
}

func (s *discoveryState) runDiscover(hook progressHook) error {
	// 发现阶段是顺序的：第 k 个模体依赖前 k-1 个模体的占用区间。
	var blocked map[int][]interval
	if s.discoverDone > 0 {
		blocked = s.blockedUntil(s.discoverDone - 1)
	}
	for k := s.discoverDone; k < s.cfg.ExpectMotifs; k++ {
		windows := enumerateWindows(s.working, s.cfg.MotifLength, blocked)
		if len(windows) < len(s.seqs) {
			break // 不足以让每条序列至少贡献一个候选位点
		}
		model, err := fitOneMotif(windows, s.cfg.MotifLength, s.cfg.Background,
			s.cfg.Pseudocount, s.cfg.ConvergenceTol, s.cfg.InitRounds,
			s.cfg.MaxEMIterations,
			s.cfg.Seed, globalStream(highDiscover, uint64(k)), "observed")
		if err != nil {
			return err
		}
		s.models[k] = model
		s.discoverDone = k + 1
		blocked = updateBlocked(blocked, model)
		if hook != nil {
			if err := hook("discover", s.discoverDone, s.cfg.ExpectMotifs); err != nil {
				return err
			}
		}
	}

	s.discoverComplete = true
	s.models = compactModels(s.models)
	s.nullScores = make([][]float64, len(s.models))
	for k := range s.models {
		s.nullScores[k] = []float64{}
	}
	// 发现阶段结束再落一次盘，记录实际模体数量与阶段完成标志；
	// 这对“窗口不足提前停止”的情况尤其重要。
	if hook != nil {
		if err := hook("discover", s.discoverDone, s.cfg.ExpectMotifs); err != nil {
			return err
		}
	}
	return nil
}

func (s *discoveryState) runShuffles(hook progressHook) error {
	for b := s.shuffleDone; b < s.cfg.Shuffles; b++ {
		shuffled := shuffleAll(s.working,
			streamRand(s.cfg.Seed, globalStream(highShuffle, uint64(b))))
		var shuffleBlocked map[int][]interval
		for k := range s.models {
			nullWindows := enumerateWindows(shuffled, s.cfg.MotifLength, shuffleBlocked)
			if len(nullWindows) < len(shuffled) {
				break
			}
			nullModel, err := fitOneMotif(nullWindows, s.cfg.MotifLength, s.cfg.Background,
				s.cfg.Pseudocount, s.cfg.ConvergenceTol, s.cfg.InitRounds,
				s.cfg.MaxEMIterations,
				s.cfg.Seed,
				globalStream(highShuffleEM, uint64(b)*uint64(s.cfg.ExpectMotifs)+uint64(k)),
				"null")
			if err != nil {
				return err
			}
			s.nullScores[k] = append(s.nullScores[k],
				motifStatistic(shuffled, nullModel.matrix))
			shuffleBlocked = updateBlocked(shuffleBlocked, nullModel)
		}
		s.shuffleDone = b + 1
		if hook != nil {
			if err := hook("shuffle", s.shuffleDone, s.cfg.Shuffles); err != nil {
				return err
			}
		}
	}
	return nil
}

func (s *discoveryState) finalize() (*Result, error) {
	if s.result != nil {
		return s.result, nil
	}
	observed := make([]float64, len(s.models))
	for k, model := range s.models {
		observed[k] = motifStatistic(s.working, model.matrix)
	}
	pvalues := make([]float64, len(s.models))
	for k := range s.models {
		pvalues[k] = empiricalPValue(observed[k], s.nullScores[k])
	}
	adjusted := benjaminiHochberg(pvalues)

	result := &Result{Motifs: make([]Motif, 0, len(s.models))}
	for k, model := range s.models {
		sites := filterSites(model, s.working, *s.cfg.ScoreThreshold)
		result.Motifs = append(result.Motifs, Motif{
			Index:          k,
			Matrix:         model.matrix,
			Sites:          sites,
			Information:    informationContent(model.matrix),
			Statistic:      observed[k],
			NullCount:      len(s.nullScores[k]),
			PValue:         pvalues[k],
			AdjustedPValue: adjusted[k],
			Significant:    adjusted[k] <= *s.cfg.SignificanceAlpha,
		})
	}
	s.result = result
	return result, nil
}

// blockedUntil 根据已发现的前 lastK+1 个模体重建占用区间集合。
func (s *discoveryState) blockedUntil(lastK int) map[int][]interval {
	var blocked map[int][]interval
	for k := 0; k <= lastK && k < len(s.models); k++ {
		if s.models[k] != nil {
			blocked = updateBlocked(blocked, s.models[k])
		}
	}
	return blocked
}

// 随机流高级编号，见 rng.go 中的两级层级说明。
const (
	highDiscover  = 1
	highShuffle   = 2
	highShuffleEM = 3
)

func compactModels(models []*emModel) []*emModel {
	out := models[:0]
	for _, m := range models {
		if m != nil {
			out = append(out, m)
		}
	}
	return out
}

func updateBlocked(blocked map[int][]interval, model *emModel) map[int][]interval {
	if blocked == nil {
		blocked = make(map[int][]interval)
	}
	for _, site := range model.sites {
		blocked[site.SeqIndex] = append(blocked[site.SeqIndex], interval{
			start: site.Offset,
			end:   site.Offset + model.matrix.Width,
		})
	}
	return blocked
}

// motifStatistic 计算模体在给定序列集上的统计量：
// 每条序列取最佳归一化得分，再对序列求和。
func motifStatistic(seqs [][]byte, matrix *PWM) float64 {
	total := 0.0
	width := matrix.Width
	for _, seq := range seqs {
		best := 0.0
		found := false
		for offset := 0; offset+width <= len(seq); offset++ {
			window := seq[offset : offset+width]
			valid := true
			for j := 0; j < width; j++ {
				if window[j] == 'N' {
					valid = false
					break
				}
			}
			if !valid {
				continue
			}
			score := matrix.NormalizedScore(window)
			if !found || score > best {
				best = score
				found = true
			}
		}
		total += best
	}
	return total
}

// filterSites 仅保留归一化得分达到阈值的位点，输出按序列与位置排序。
func filterSites(model *emModel, seqs [][]byte, threshold float64) []Site {
	width := model.matrix.Width
	var sites []Site
	for si, seq := range seqs {
		for offset := 0; offset+width <= len(seq); offset++ {
			window := seq[offset : offset+width]
			valid := true
			for j := 0; j < width; j++ {
				if window[j] == 'N' {
					valid = false
					break
				}
			}
			if !valid {
				continue
			}
			score := model.matrix.NormalizedScore(window)
			if score >= threshold {
				sites = append(sites, Site{SeqIndex: si, Offset: offset, Score: score})
			}
		}
	}
	return sites
}

// informationContent 返回 PWM 的平均信息含量（比特/位点）。
func informationContent(m *PWM) float64 {
	total := 0.0
	for _, row := range m.LogOdds {
		for b := 0; b < 4; b++ {
			prior := m.bg.Freq[b]
			if prior <= 0 {
				prior = 0.25
			}
			prob := prior * math.Exp(row[b])
			if prob > 0 {
				total += prob * math.Log2(prob/prior)
			}
		}
	}
	return total / float64(m.Width)
}
