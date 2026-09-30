package motif

import (
	"fmt"
	"math"
)

// 核苷酸索引约定：A=0, C=1, G=2, T=3。

var baseIndex [256]int8

func init() {
	for i := range baseIndex {
		baseIndex[i] = -1
	}
	baseIndex['A'] = 0
	baseIndex['C'] = 1
	baseIndex['G'] = 2
	baseIndex['T'] = 3
}

// UniformBackground 返回均匀背景 (0.25,0.25,0.25,0.25)。
func UniformBackground() Background {
	return Background{0.25, 0.25, 0.25, 0.25}
}

// Validate 校验背景概率是否非负且总和为 1。
func (b Background) Validate() error {
	sum := 0.0
	for _, p := range b {
		if p < 0 || math.IsNaN(p) || math.IsInf(p, 0) {
			return fmt.Errorf("motif: 背景概率非法：%v", b)
		}
		sum += p
	}
	if math.Abs(sum-1) > 1e-9 {
		return fmt.Errorf("motif: 背景概率总和必须为 1，实际为 %g", sum)
	}
	return nil
}

// Background 是四核苷酸背景概率分布，顺序 A,C,G,T，总和为 1。
type Background [4]float64

// PWM 是长度为 L 的位置权重矩阵。
// Weight[i][b] 为位置 i 上碱基 b 的加性对数似然比权重。
type PWM struct {
	Length     int
	Weight     [][4]float64
	Background Background
	Consensus  string
	Info       []float64
}

// PWMConfig 控制 PWM 构建与打分参数。
type PWMConfig struct {
	Pseudocount float64
	Threshold   float64
}

// Hit 描述一次矩阵命中。
type Hit struct {
	SeqIndex int
	Start    int
	Score    float64
}

// DefaultPWMConfig 返回默认的矩阵构建与打分配置。
func DefaultPWMConfig() PWMConfig {
	return PWMConfig{
		Pseudocount: 0.8,
		Threshold:   0.8,
	}
}

// BuildPWM 依据对齐片段与背景构建带伪计数的位置权重矩阵。
func BuildPWM(sites [][]byte, bg Background, cfg PWMConfig) (*PWM, error) {
	if len(sites) == 0 {
		return nil, fmt.Errorf("motif: 至少需要一个对齐片段才能构建 PWM")
	}
	if err := bg.Validate(); err != nil {
		return nil, err
	}
	length := len(sites[0])
	if length == 0 {
		return nil, fmt.Errorf("motif: 片段长度必须为正数")
	}
	pseudo := cfg.Pseudocount
	if pseudo < 0 {
		return nil, fmt.Errorf("motif: 伪计数不能为负数")
	}
	counts := make([][4]float64, length)
	for si, site := range sites {
		if len(site) != length {
			return nil, fmt.Errorf("motif: 第 %d 个片段长度 %d 与期望长度 %d 不一致", si, len(site), length)
		}
		for i, c := range site {
			b := baseIndex[c]
			if b < 0 {
				return nil, fmt.Errorf("motif: 片段 %d 位置 %d 含非法碱基 %q", si, i, c)
			}
			counts[i][b]++
		}
	}

	n := float64(len(sites))
	weight := make([][4]float64, length)
	info := make([]float64, length)
	cons := make([]byte, length)
	for i := 0; i < length; i++ {
		bestBase, bestProb := 0, -1.0
		for b := 0; b < 4; b++ {
			freq := (counts[i][b] + pseudo*bg[b]) / (n + pseudo)
			weight[i][b] = math.Log(freq / bg[b])
			info[i] += freq * (math.Log2(freq) + 2)
			if freq > bestProb {
				bestProb = freq
				bestBase = b
			}
		}
		cons[i] = "ACGT"[bestBase]
	}

	m := &PWM{
		Length:     length,
		Weight:     weight,
		Background: bg,
		Consensus:  string(cons),
		Info:       info,
	}
	return m, nil
}

// BuildWeightedPWM 依据带权重的候选位点构建 PWM。
// data[i] 为第 i 条序列；starts[i][p] 为第 p 个候选起点；weights[i][p]
// 为该候选起点的后验责任（E 步输出），三者长度关系为
// len(weights[i]) == len(starts[i])。
func BuildWeightedPWM(length int, data [][]byte, starts [][]int, weights [][]float64, bg Background, cfg PWMConfig) (*PWM, error) {
	if len(data) == 0 || len(data) != len(weights) || len(data) != len(starts) {
		return nil, fmt.Errorf("motif: 加权 PWM 输入不合法")
	}
	if err := bg.Validate(); err != nil {
		return nil, err
	}
	pseudo := cfg.Pseudocount
	if pseudo < 0 {
		return nil, fmt.Errorf("motif: 伪计数不能为负数")
	}
	if length <= 0 {
		return nil, fmt.Errorf("motif: 加权 PWM 长度必须为正数")
	}

	counts := make([][4]float64, length)
	totalW := 0.0
	for i, seq := range data {
		if len(weights[i]) == 0 {
			continue
		}
		for p, w := range weights[i] {
			if w <= 0 {
				continue
			}
			origin := starts[i][p]
			for pos := 0; pos < length; pos++ {
				b := baseIndex[seq[origin+pos]]
				if b < 0 {
					return nil, fmt.Errorf("motif: 序列 %d 位置 %d 含非法碱基", i, p+pos)
				}
				counts[pos][b] += w
			}
			totalW += w
		}
	}
	if totalW <= 0 {
		return nil, fmt.Errorf("motif: 后验权重总和为零")
	}

	weight := make([][4]float64, length)
	info := make([]float64, length)
	cons := make([]byte, length)
	for i := 0; i < length; i++ {
		bestBase, bestProb := 0, -1.0
		for b := 0; b < 4; b++ {
			freq := (counts[i][b] + pseudo*bg[b]) / (totalW + pseudo)
			weight[i][b] = math.Log(freq / bg[b])
			info[i] += freq * (math.Log2(freq) + 2)
			if freq > bestProb {
				bestProb = freq
				bestBase = b
			}
		}
		cons[i] = "ACGT"[bestBase]
	}
	return &PWM{
		Length:     length,
		Weight:     weight,
		Background: bg,
		Consensus:  string(cons),
		Info:       info,
	}, nil
}

// Score 返回单个长度为 L 的片段相对背景模型的归一化得分。
func (m *PWM) Score(site []byte) float64 {
	s := 0.0
	for i, c := range site {
		b := baseIndex[c]
		if b < 0 {
			return math.Inf(-1)
		}
		s += m.Weight[i][b]
	}
	return s
}

// MaxScore 返回矩阵的理论最大得分。
func (m *PWM) MaxScore() float64 {
	return m.extreme(1)
}

// MinScore 返回矩阵的理论最小得分。
func (m *PWM) MinScore() float64 {
	return m.extreme(-1)
}

// NormalizedScore 返回归一化到 [0,1] 区间的得分。
func (m *PWM) NormalizedScore(site []byte) float64 {
	raw := m.Score(site)
	max := m.MaxScore()
	min := m.MinScore()
	if max <= min {
		return 0
	}
	v := (raw - min) / (max - min)
	if v < 0 {
		return 0
	}
	if v > 1 {
		return 1
	}
	return v
}

// Scan 在一条序列中滑动扫描，返回得分达到阈值的所有命中。
// threshold 为归一化阈值（[0,1]，与 NormalizedScore 比较）；
// 传入小于 0 的值可获取每个位置的原始最大命中。
func (m *PWM) Scan(seq []byte, threshold float64) []Hit {
	var hits []Hit
	if len(seq) < m.Length {
		return hits
	}
	max := m.MaxScore()
	min := m.MinScore()
	denom := max - min
	for start := 0; start+m.Length <= len(seq); start++ {
		raw := m.Score(seq[start : start+m.Length])
		if math.IsInf(raw, -1) {
			continue
		}
		norm := 0.0
		if denom > 0 {
			norm = (raw - min) / denom
		}
		if norm >= threshold {
			hits = append(hits, Hit{Start: start, Score: norm})
		}
	}
	return hits
}

// ScanAll 在所有序列上扫描矩阵命中。
func (m *PWM) ScanAll(seqs []Seq, threshold float64) []Hit {
	var hits []Hit
	for i, s := range seqs {
		locals := m.Scan(s.Data, threshold)
		for j := range locals {
			locals[j].SeqIndex = i
		}
		hits = append(hits, locals...)
	}
	return hits
}

func (m *PWM) extreme(sign float64) float64 {
	s := 0.0
	for i := 0; i < m.Length; i++ {
		v := m.Weight[i][0]
		for _, w := range m.Weight[i][1:] {
			if sign > 0 && w > v || sign < 0 && w < v {
				v = w
			}
		}
		s += v
	}
	return s
}
