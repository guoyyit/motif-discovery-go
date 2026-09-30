package motif

import "fmt"

// Config 是模体发现流程的完整参数配置。
type Config struct {
	// MotifLength 为待发现模体长度 L。
	MotifLength int
	// ExpectedCount 为期望输出的模体数量 K。
	ExpectedCount int
	// InitRounds 为多轮随机初始化的轮数。
	InitRounds int
	// ConvergenceTol 为 EM 收敛容差（对数似然相对变化）。
	ConvergenceTol float64
	// MaxEMIterations 为单轮 EM 最大迭代次数。
	MaxEMIterations int
	// PWM 控制矩阵构建与打分阈值。
	PWM PWMConfig
	// Filter 控制低复杂度/重复区域处理。
	Filter FilterConfig
	// Significance 控制随机重排显著性检验。
	Significance SigConfig
	// Seed 为控制所有随机过程的显式种子。
	Seed uint64
}

// DefaultConfig 返回合理的默认配置。
func DefaultConfig() Config {
	return Config{
		MotifLength:     10,
		ExpectedCount:   1,
		InitRounds:      20,
		ConvergenceTol:  1e-4,
		MaxEMIterations: 100,
		PWM:             DefaultPWMConfig(),
		Filter:          DefaultFilterConfig(),
		Significance:    DefaultSigConfig(),
		Seed:            1,
	}
}

// Validate 校验配置取值，非法时返回错误。
func (c Config) Validate() error {
	if c.MotifLength <= 0 {
		return fmt.Errorf("motif: MotifLength 必须为正数，实际为 %d", c.MotifLength)
	}
	if c.ExpectedCount <= 0 {
		return fmt.Errorf("motif: ExpectedCount 必须为正数，实际为 %d", c.ExpectedCount)
	}
	if c.InitRounds <= 0 {
		return fmt.Errorf("motif: InitRounds 必须为正数，实际为 %d", c.InitRounds)
	}
	if c.ConvergenceTol <= 0 {
		return fmt.Errorf("motif: ConvergenceTol 必须为正数，实际为 %g", c.ConvergenceTol)
	}
	if c.MaxEMIterations <= 0 {
		return fmt.Errorf("motif: MaxEMIterations 必须为正数，实际为 %d", c.MaxEMIterations)
	}
	if c.PWM.Pseudocount < 0 {
		return fmt.Errorf("motif: PWM.Pseudocount 不能为负数")
	}
	if c.PWM.Threshold < 0 || c.PWM.Threshold > 1 {
		return fmt.Errorf("motif: PWM.Threshold 必须落在 [0,1]，实际为 %g", c.PWM.Threshold)
	}
	if c.Significance.Permutations < 0 {
		return fmt.Errorf("motif: Significance.Permutations 不能为负数")
	}
	if c.Significance.Alpha <= 0 || c.Significance.Alpha > 1 {
		return fmt.Errorf("motif: Significance.Alpha 必须落在 (0,1]，实际为 %g", c.Significance.Alpha)
	}
	if c.Significance.Method != "BH" {
		return fmt.Errorf("motif: Significance.Method 仅支持 \"BH\"，实际为 %q", c.Significance.Method)
	}
	return nil
}

// Motif 是一个发现得到的模体结果。
type Motif struct {
	Index     int
	Matrix    *PWM
	Consensus string
	Score     float64
	// Enrichment 为统计检验统计量：所有序列中归一化得分达到 PWM.Threshold
	// 的位点比例（每条序列最多计一次），直接衡量模体富集程度。
	Enrichment  float64
	Sites       []SiteOccurrence
	Information float64
	PValue      float64
	AdjustedP   float64
	Significant bool
}

// SiteOccurrence 描述模体在某条序列上的一个发生位点。
type SiteOccurrence struct {
	SeqIndex int
	Start    int
	Score    float64
}

// Result 是一次发现任务的完整输出。
type Result struct {
	Seed        uint64
	Motifs      []Motif
	Background  Background
	SequenceCnt int
	TotalLen    int
}
