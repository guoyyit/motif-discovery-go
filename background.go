package motif

import (
	"fmt"
	"math"
)

// Background 是 A/C/G/T 的一阶背景频率模型。
type Background struct {
	Freq [4]float64
}

// baseIndex 返回 A/C/G/T 的下标 0..3，其余字符（含 N）返回 -1。
func baseIndex(base byte) int {
	switch base {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	default:
		return -1
	}
}

// DefaultBackground 返回均匀背景模型。
func DefaultBackground() Background {
	return Background{Freq: [4]float64{0.25, 0.25, 0.25, 0.25}}
}

// EstimateBackground 从序列集合估计背景频率。
// 未知碱基 N 不计入统计；若没有任何已知碱基则退化为均匀模型。
func EstimateBackground(seqs []Sequence) Background {
	var counts [4]float64
	total := 0.0
	for _, seq := range seqs {
		for _, base := range seq.Data {
			if idx := baseIndex(base); idx >= 0 {
				counts[idx]++
				total++
			}
		}
	}
	if total == 0 {
		return DefaultBackground()
	}
	bg := Background{}
	for i := range counts {
		bg.Freq[i] = counts[i] / total
	}
	return bg
}

// Prob 返回碱基在背景模型下的概率，未知碱基返回 0。
func (b Background) Prob(base byte) float64 {
	if idx := baseIndex(base); idx >= 0 {
		return b.Freq[idx]
	}
	return 0
}

// Validate 检查背景模型是否为合法概率分布。
func (b Background) Validate() error {
	sum := 0.0
	for i, p := range b.Freq {
		if math.IsNaN(p) || p < 0 {
			return fmt.Errorf("background: frequency for base %d is negative or NaN: %v", i, p)
		}
		sum += p
	}
	if math.Abs(sum-1) > 1e-9 {
		return fmt.Errorf("background: frequencies must sum to 1, got %v (sum=%g)", b.Freq, sum)
	}
	return nil
}
