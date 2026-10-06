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
	// 观测拟合与零分布拟合使用完全相同的轮数。
	InitRounds int
	// ConvergenceTol 是 EM 对数似然的收敛容差（绝对变化量）。
	ConvergenceTol float64
	// MaxIterations 是单次 EM 拟合的最大迭代步数；达到上限仍未收敛时
	// 该轮拟合以不收敛处理（见 ErrNotConverged）。
	MaxIterations int
	// ScoreThreshold 是位点归一化得分阈值（(0,1]），低于该值的位点不计入。
	ScoreThreshold float64
	// SignificanceAlpha 是经多重检验校正后的显著性阈值（FDR，(0,1]）。
	SignificanceAlpha float64
	// Shuffles 是显著性零分布的重排次数 B，最小经验 p 值为 1/(B+1)。
	Shuffles int
	// Pseudocount 是 PWM 的伪计数强度。
	Pseudocount float64
	// Background 为零值（全零）时自动从输入序列估计。
	Background Background
	// MaskLowComplexity 为 true 时先掩码低复杂度区域，简单重复不会成为模体。
	// 掩码自始即为默认行为，该字段保留用于兼容显式开启的旧配置。
	MaskLowComplexity bool
	// NoMaskLowComplexity 显式关闭低复杂度掩码防护。默认（未设置）时掩码
	// 始终开启；只有明确设置该字段才恢复不掩码的旧行为。
	// 与 MaskLowComplexity 同时设置属于配置冲突，校验会报错。
	NoMaskLowComplexity bool
}

// Site 是模体在某条序列中的一个出现位置。
type Site struct {
	SeqIndex int
	Offset   int
	Score    float64
}

// Motif 是单个发现模体及其统计结论。
type Motif struct {
	Index          int
	Matrix         *PWM
	Sites          []Site
	Information    float64
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
	defaultMaxIterations  = 500
	defaultScoreThreshold = 0.7
	defaultAlpha          = 0.05
	defaultShuffles       = 200
	defaultPseudocount    = 0.5
)

// Validate 检查配置合法性并填充默认值，返回填充后的配置。
//
// 零值字段视为“未配置”，自动填入默认值；显式给出的越界值
// （负数计数、超出合法区间的阈值、NaN 等）直接报错，说明字段名与
// 合法取值范围，绝不静默替换为默认值后继续运行。
func (c Config) Validate() (Config, error) {
	// 计数参数：零值填默认，负值报错。每个字段独立判断。
	if c.MotifLength < 0 {
		return c, fmt.Errorf("config: MotifLength must be >= 4 (0 selects the default %d), got %d",
			defaultMotifLength, c.MotifLength)
	}
	if c.MotifLength == 0 {
		c.MotifLength = defaultMotifLength
	}
	if c.ExpectMotifs < 0 {
		return c, fmt.Errorf("config: ExpectMotifs must be >= 1 (0 selects the default %d), got %d",
			defaultExpectMotifs, c.ExpectMotifs)
	}
	if c.ExpectMotifs == 0 {
		c.ExpectMotifs = defaultExpectMotifs
	}
	if c.InitRounds < 0 {
		return c, fmt.Errorf("config: InitRounds must be >= 1 (0 selects the default %d), got %d",
			defaultInitRounds, c.InitRounds)
	}
	if c.InitRounds == 0 {
		c.InitRounds = defaultInitRounds
	}
	if c.MaxIterations < 0 {
		return c, fmt.Errorf("config: MaxIterations must be >= 1 (0 selects the default %d), got %d",
			defaultMaxIterations, c.MaxIterations)
	}
	if c.MaxIterations == 0 {
		c.MaxIterations = defaultMaxIterations
	}
	if c.Shuffles < 0 {
		return c, fmt.Errorf("config: Shuffles must be >= 1 (0 selects the default %d), got %d",
			defaultShuffles, c.Shuffles)
	}
	if c.Shuffles == 0 {
		c.Shuffles = defaultShuffles
	}

	// 浮点参数：零值填默认，负值、越界值与 NaN 报错。
	if math.IsNaN(c.ConvergenceTol) || c.ConvergenceTol < 0 {
		return c, fmt.Errorf("config: ConvergenceTol must be > 0 (0 selects the default %g), got %g",
			defaultConvergenceTol, c.ConvergenceTol)
	}
	if c.ConvergenceTol == 0 {
		c.ConvergenceTol = defaultConvergenceTol
	}
	if math.IsNaN(c.ScoreThreshold) || c.ScoreThreshold < 0 || c.ScoreThreshold > 1 {
		return c, fmt.Errorf("config: ScoreThreshold must be within (0,1] (0 selects the default %g), got %g",
			defaultScoreThreshold, c.ScoreThreshold)
	}
	if c.ScoreThreshold == 0 {
		c.ScoreThreshold = defaultScoreThreshold
	}
	if math.IsNaN(c.SignificanceAlpha) || c.SignificanceAlpha < 0 || c.SignificanceAlpha > 1 {
		return c, fmt.Errorf("config: SignificanceAlpha must be within (0,1] (0 selects the default %g), got %g",
			defaultAlpha, c.SignificanceAlpha)
	}
	if c.SignificanceAlpha == 0 {
		c.SignificanceAlpha = defaultAlpha
	}
	if math.IsNaN(c.Pseudocount) || c.Pseudocount < 0 {
		return c, fmt.Errorf("config: Pseudocount must be >= 0 (0 selects the default %g), got %g",
			defaultPseudocount, c.Pseudocount)
	}
	if c.Pseudocount == 0 {
		c.Pseudocount = defaultPseudocount
	}

	// 填充默认值后再做跨字段与下界校验。
	if c.MotifLength < 4 {
		return c, fmt.Errorf("config: MotifLength must be >= 4, got %d", c.MotifLength)
	}
	if c.MaskLowComplexity && c.NoMaskLowComplexity {
		return c, fmt.Errorf("config: MaskLowComplexity and NoMaskLowComplexity are conflicting; set at most one")
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

// maskEnabled 报告低复杂度掩码是否生效：默认开启，只有显式设置
// NoMaskLowComplexity 才关闭；MaskLowComplexity 为历史兼容字段。
func maskEnabled(cfg Config) bool {
	return !cfg.NoMaskLowComplexity
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
		if maskEnabled(cfg) {
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
			s.cfg.MaxIterations, s.cfg.Seed, globalStream(highDiscover, uint64(k)))
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
				s.cfg.MaxIterations, s.cfg.Seed,
				globalStream(highShuffleEM, uint64(b)*uint64(s.cfg.ExpectMotifs)+uint64(k)))
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
		sites := filterSites(model, s.working, s.cfg.ScoreThreshold)
		result.Motifs = append(result.Motifs, Motif{
			Index:          k,
			Matrix:         model.matrix,
			Sites:          sites,
			Information:    informationContent(model.matrix),
			PValue:         pvalues[k],
			AdjustedPValue: adjusted[k],
			Significant:    adjusted[k] <= s.cfg.SignificanceAlpha,
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
