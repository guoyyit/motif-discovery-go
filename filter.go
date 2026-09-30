package motif

import "math"

// FilterConfig 控制低复杂度区域与简单重复序列的处理。
type FilterConfig struct {
	// Enabled 为 false 时不过滤任何区域。
	Enabled bool
	// WindowSize 为滑窗熵检测的窗口长度。
	WindowSize int
	// MaxEntropy 为窗口允许的最大香农熵（bit），低于该值标记为低复杂度。
	MaxEntropy float64
	// MonomerRunLength 达到该长度的单碱基连续延伸被标记。
	MonomerRunLength int
	// DinucRepeatLen 达到该长度（bp）的双碱基串联重复被标记。
	DinucRepeatLen int
}

// DefaultFilterConfig 返回默认的低复杂度过滤配置。
func DefaultFilterConfig() FilterConfig {
	return FilterConfig{
		Enabled:          true,
		WindowSize:       12,
		MaxEntropy:       1.0,
		MonomerRunLength: 8,
		DinucRepeatLen:   10,
	}
}

// MaskedRegion 描述一段被标记为低复杂度/重复的区域。
type MaskedRegion struct {
	Start int
	End   int
	Kind  string
}

// ComplexityFilter 对序列标记低复杂度与简单重复区域。
type ComplexityFilter struct {
	cfg FilterConfig
}

// NewComplexityFilter 构造低复杂度过滤器。
func NewComplexityFilter(cfg FilterConfig) *ComplexityFilter {
	if cfg.Enabled && cfg.WindowSize <= 0 {
		d := DefaultFilterConfig()
		cfg.WindowSize = d.WindowSize
		cfg.MaxEntropy = d.MaxEntropy
		cfg.MonomerRunLength = d.MonomerRunLength
		cfg.DinucRepeatLen = d.DinucRepeatLen
	}
	return &ComplexityFilter{cfg: cfg}
}

// MaskRegions 返回序列中的低复杂度与重复区域（半开区间）。
func (f *ComplexityFilter) MaskRegions(seq []byte) []MaskedRegion {
	if !f.cfg.Enabled || len(seq) == 0 {
		return nil
	}
	n := len(seq)
	mask := make([]bool, n)
	kind := make([]string, n)

	mark := func(start, end int, label string) {
		for i := start; i < end; i++ {
			if !mask[i] {
				mask[i] = true
				kind[i] = label
			}
		}
	}

	// 单碱基连续延伸，例如 AAAAAAAA。
	if f.cfg.MonomerRunLength > 0 {
		i := 0
		for i < n {
			j := i + 1
			for j < n && seq[j] == seq[i] {
				j++
			}
			if j-i >= f.cfg.MonomerRunLength {
				mark(i, j, "monomer-run")
			}
			i = j
		}
	}

	// 周期 2..4 的串联重复，例如 ATATATAT、ACGACGACG、ACGTACGT。
	if f.cfg.DinucRepeatLen > 0 {
		for period := 2; period <= 4 && period <= n; period++ {
			flushTandem := func(start, end int) {
				repeats := end - start
				if repeats < f.cfg.DinucRepeatLen || repeats < 3*period {
					return
				}
				if !unitNonDegenerate(seq[start : start+period]) {
					return
				}
				for k := start + period; k < end; k++ {
					if seq[k] != seq[start+(k-start)%period] {
						return
					}
				}
				for u := start; u+period <= end; u += period {
					if !unitNonDegenerate(seq[u : u+period]) {
						return
					}
				}
				if period > 2 && explainedBySmallerPeriod(seq, start, end, period) {
					return
				}
				mark(start, end, "tandem-repeat")
			}
			runStart := -1
			for i := period; i <= n; i++ {
				match := i < n && seq[i] == seq[i-period]
				if match {
					if runStart < 0 {
						// 至少需要连续两个匹配位置才确认游程起点，
						// 排除 ACGT|AAAA 边界处单个巧合匹配。
						if i+1 < n && seq[i+1] == seq[i+1-period] {
							runStart = i
						}
					}
					continue
				}
				if runStart >= 0 {
					flushTandem(runStart-period, i)
					runStart = -1
				}
			}
		}
	}

	// 滑窗香农熵。
	if f.cfg.WindowSize > 0 && f.cfg.MaxEntropy > 0 {
		w := f.cfg.WindowSize
		for start := 0; start+w <= n; start++ {
			var counts [4]int
			ok := true
			for k := start; k < start+w; k++ {
				b := baseIndex[seq[k]]
				if b < 0 {
					ok = false
					break
				}
				counts[b]++
			}
			if ok && windowEntropy(counts[:], w) < f.cfg.MaxEntropy {
				mark(start, start+w, "low-entropy")
			}
		}
	}

	return maskToRegions(mask, kind)
}

// IsMasked 判断位置 pos 是否落在被标记区域内。
func (f *ComplexityFilter) IsMasked(regions []MaskedRegion, pos int) bool {
	for _, r := range regions {
		if pos >= r.Start && pos < r.End {
			return true
		}
	}
	return false
}

// AllowedStarts 返回给定序列中允许放置模体起始位点的位置集合。
func (f *ComplexityFilter) AllowedStarts(seq []byte, motifLen int) []int {
	if motifLen <= 0 || len(seq) < motifLen {
		return nil
	}
	regions := f.MaskRegions(seq)
	out := make([]int, 0, len(seq)-motifLen+1)
	for start := 0; start+motifLen <= len(seq); start++ {
		end := start + motifLen
		blocked := false
		for _, r := range regions {
			if r.End <= start {
				continue
			}
			if r.Start >= end {
				break
			}
			blocked = true
			break
		}
		if !blocked {
			out = append(out, start)
		}
	}
	return out
}

// SitePasses 判断候选片段本身是否通过简单重复检查。
func (f *ComplexityFilter) SitePasses(site []byte) bool {
	if !f.cfg.Enabled || len(site) == 0 {
		return true
	}
	n := len(site)
	var counts [4]int
	for _, c := range site {
		b := baseIndex[c]
		if b < 0 {
			return false
		}
		counts[b]++
	}
	// 单碱基占比几乎为 1。
	if n >= 5 {
		for _, c := range counts {
			if c == n {
				return false
			}
		}
	}
	// 周期 1..4 的简单重复：至少覆盖 2 个单元且 80% 以上位点符合周期。
	for period := 1; period <= 4 && 2*period <= n; period++ {
		matches := 0
		total := 0
		for i := period; i < n; i++ {
			total++
			if site[i] == site[i-period] {
				matches++
			}
		}
		if total > 0 && float64(matches)/float64(total) >= 0.8 {
			return false
		}
	}
	// 片段自身低熵。
	if f.cfg.MaxEntropy > 0 && windowEntropy(counts[:], n) < f.cfg.MaxEntropy {
		return false
	}
	return true
}

func windowEntropy(counts []int, total int) float64 {
	h := 0.0
	for _, c := range counts {
		if c == 0 {
			continue
		}
		p := float64(c) / float64(total)
		h -= p * math.Log2(p)
	}
	return h
}

// unitNonDegenerate 判断周期单元内至少含两种碱基，
// 排除 AAAAAA 被当作任意周期重复之类的退化情况。
func unitNonDegenerate(unit []byte) bool {
	for d := 1; d < len(unit); d++ {
		if unit[d] != unit[0] {
			return true
		}
	}
	return false
}

// explainedBySmallerPeriod 判断区间是否已由 1..period-1 的更短周期完全解释，
// 避免 ACACACAC 同时以周期 2 和周期 4 被报告。
func explainedBySmallerPeriod(seq []byte, start, end, period int) bool {
	for p := 1; p < period; p++ {
		ok := true
		for i := start + p; i < end; i++ {
			if seq[i] != seq[i-p] {
				ok = false
				break
			}
		}
		if ok {
			return true
		}
	}
	return false
}

func maskToRegions(mask []bool, kind []string) []MaskedRegion {
	var regions []MaskedRegion
	for i := 0; i < len(mask); i++ {
		if !mask[i] {
			continue
		}
		start := i
		label := kind[i]
		for i < len(mask) && mask[i] && kind[i] == label {
			i++
		}
		regions = append(regions, MaskedRegion{Start: start, End: i, Kind: label})
		i--
	}
	return regions
}
