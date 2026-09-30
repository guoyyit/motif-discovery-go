package motif

import (
	"math"
	"sort"
)

// empiricalPValue 根据重排空分布计算单侧经验 p 值：
// p = (1 + #{null >= observed}) / (1 + len(nullScores))。
// +1 修正（Davison & Hinkley）保证 p 值严格大于 0，有限次重排下
// 可达到的最小 p 值为 1/(B+1)。
func empiricalPValue(observed float64, nullScores []float64) float64 {
	b := len(nullScores)
	exceeds := 0
	for _, null := range nullScores {
		if null >= observed {
			exceeds++
		}
	}
	return float64(1+exceeds) / float64(1+b)
}

// benjaminiHochberg 执行 Benjamini-Hochberg FDR 多重检验校正。
// 输入是 m 个原始 p 值（对应同一发现任务中的全部候选模体），
// 返回等长的校正后 q 值；输出顺序与输入一致。
// NaN p 值按 1.0 处理。
func benjaminiHochberg(pvalues []float64) []float64 {
	m := len(pvalues)
	adjusted := make([]float64, m)
	if m == 0 {
		return adjusted
	}

	type entry struct {
		index int
		p     float64
	}
	entries := make([]entry, m)
	for i, p := range pvalues {
		if math.IsNaN(p) {
			p = 1
		}
		entries[i] = entry{index: i, p: p}
	}
	sort.Slice(entries, func(i, j int) bool {
		if entries[i].p != entries[j].p {
			return entries[i].p < entries[j].p
		}
		return entries[i].index < entries[j].index
	})

	// 按秩从大到小做单调性保证：q_(k) = min over j>=k of m/j * p_(j)。
	prev := 1.0
	for rank := m; rank >= 1; rank-- {
		e := entries[rank-1]
		q := float64(m) / float64(rank) * e.p
		if q > prev {
			q = prev
		}
		if q > 1 {
			q = 1
		}
		adjusted[e.index] = q
		prev = q
	}
	return adjusted
}
