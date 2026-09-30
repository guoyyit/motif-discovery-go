package motif

import (
	"fmt"
	"math"
)

// PWM 是以对数几率表示的位置权重矩阵，行是模体位置，列是 A/C/G/T。
type PWM struct {
	Width   int
	LogOdds [][]float64
	bg      Background

	min float64
	max float64
}

// BuildPWM 从一组等长比对片段构建 PWM。
// 每个位点先统计碱基计数并加入 pseudocount*背景频率 的伪计数，
// 再取相对背景模型的自然对数几率。含 N 的片段会被拒绝。
func BuildPWM(aligned [][]byte, bg Background, pseudocount float64) (*PWM, error) {
	if len(aligned) == 0 {
		return nil, fmt.Errorf("pwm: no aligned windows provided")
	}
	width := len(aligned[0])
	if width == 0 {
		return nil, fmt.Errorf("pwm: aligned windows must not be empty")
	}
	for i, window := range aligned {
		if len(window) != width {
			return nil, fmt.Errorf("pwm: window %d has length %d, want %d", i, len(window), width)
		}
		for _, base := range window {
			if baseIndex(base) < 0 {
				return nil, fmt.Errorf("pwm: window %d contains unknown base %q", i, string(base))
			}
		}
	}
	if pseudocount < 0 {
		return nil, fmt.Errorf("pwm: pseudocount must be non-negative, got %g", pseudocount)
	}

	logOdds := make([][]float64, width)
	n := float64(len(aligned))
	for pos := 0; pos < width; pos++ {
		var counts [4]float64
		for _, window := range aligned {
			counts[baseIndex(window[pos])]++
		}
		row := make([]float64, 4)
		for b := 0; b < 4; b++ {
			prior := bg.Freq[b]
			if prior <= 0 {
				prior = 0.25
			}
			freq := (counts[b] + pseudocount*prior) / (n + pseudocount)
			row[b] = math.Log(freq / prior)
		}
		logOdds[pos] = row
	}

	m := &PWM{Width: width, LogOdds: logOdds, bg: bg}
	for _, row := range logOdds {
		lo, hi := row[0], row[0]
		for _, v := range row[1:] {
			if v < lo {
				lo = v
			}
			if v > hi {
				hi = v
			}
		}
		m.min += lo
		m.max += hi
	}
	return m, nil
}

// Score 返回窗口的原始对数几率得分（各位点之和）。
// 窗口长度必须等于矩阵宽度，且不能含未知碱基 N。
func (m *PWM) Score(window []byte) float64 {
	if len(window) != m.Width {
		panic(fmt.Sprintf("pwm: window length %d does not match motif width %d", len(window), m.Width))
	}
	score := 0.0
	for pos := 0; pos < m.Width; pos++ {
		idx := baseIndex(window[pos])
		if idx < 0 {
			panic(fmt.Sprintf("pwm: cannot score unknown base %q", string(window[pos])))
		}
		score += m.LogOdds[pos][idx]
	}
	return score
}

// NormalizedScore 返回归一化到 [0,1] 的得分：
// (raw-min)/(max-min)。归一化使不同背景/长度的矩阵阈值可比较。
func (m *PWM) NormalizedScore(window []byte) float64 {
	raw := m.Score(window)
	if m.max == m.min {
		return 0
	}
	normalized := (raw - m.min) / (m.max - m.min)
	if normalized < 0 {
		return 0
	}
	if normalized > 1 {
		return 1
	}
	return normalized
}

// MaxScore 返回矩阵的最大可能原始得分。
func (m *PWM) MaxScore() float64 { return m.max }

// MinScore 返回矩阵的最小可能原始得分。
func (m *PWM) MinScore() float64 { return m.min }

// Threshold 将相对阈值（[0,1]，基于归一化得分）换算为原始得分阈值。
func (m *PWM) Threshold(relative float64) float64 {
	if relative < 0 {
		relative = 0
	}
	if relative > 1 {
		relative = 1
	}
	return m.min + relative*(m.max-m.min)
}

// Consensus 返回每一位点取最大似然碱基得到的一致性序列。
func (m *PWM) Consensus() string {
	bases := []byte{'A', 'C', 'G', 'T'}
	out := make([]byte, m.Width)
	for pos, row := range m.LogOdds {
		best := 0
		for b := 1; b < 4; b++ {
			if row[b] > row[best] {
				best = b
			}
		}
		out[pos] = bases[best]
	}
	return string(out)
}
