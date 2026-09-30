package motif

import "math"

// 低复杂度检测参数。
const (
	// maskWindow 是低复杂度检测的滑动窗口长度。
	maskWindow = 8
	// maxRepeatPeriod 检测的简单串联重复的最大周期（如周期 2 的 ATATAT）。
	maxRepeatPeriod = 4
	// homopolymerRun 判定均聚物所需的连续相同碱基长度。
	homopolymerRun = 6
	// minEntropyBits 窗口香农熵低于该值（单位 bit）时判为低复杂度。
	minEntropyBits = 1.0
)

// MaskLowComplexity 返回将低复杂度区域替换为 N 的序列副本。
// 判定为低复杂度的固定长度窗口（maskWindow）中每个碱基都会被掩码；
// 长度不足一个窗口的尾部若为均聚物同样会被掩码。
func MaskLowComplexity(seq []byte) []byte {
	masked := append([]byte(nil), seq...)
	maskHomopolymers(masked)
	n := len(masked)
	for start := 0; start+maskWindow <= n; start++ {
		window := masked[start : start+maskWindow]
		if !IsLowComplexity(window) {
			continue
		}
		lo, hi := start, start+maskWindow
		if period, segStart, ok := repeatMatch(window); ok {
			// 重复单元在整条序列中的“相位原点”：窗口起点 + 重复起点。
			origin := start + segStart
			for lo > 0 && masked[lo-1] != 'N' &&
				masked[lo-1] == masked[origin+((lo-1-origin)%period+period)%period] {
				lo--
			}
			for hi < n && masked[hi] != 'N' &&
				masked[hi] == masked[origin+((hi-origin)%period+period)%period] {
				hi++
			}
		}
		for j := lo; j < hi; j++ {
			masked[j] = 'N'
		}
	}
	return masked
}

// IsLowComplexity 判断给定窗口是否为低复杂度/简单重复片段。
// 命中以下任一条件即判为低复杂度：
//   - 存在长度达到 homopolymerRun 的均聚物（如 AAAAAA）；
//   - 存在周期 1..maxRepeatPeriod 的串联重复，且至少重复 3 次
//     （如 ATA、ATAT）；
//   - 窗口（忽略已有的 N）香农熵低于 minEntropyBits。
func IsLowComplexity(window []byte) bool {
	if len(window) == 0 {
		return false
	}
	known := 0
	for _, base := range window {
		if baseIndex(base) >= 0 {
			known++
		}
	}
	// 已被掩码（N）占多数的窗口按复杂区域处理，避免掩码向相邻正常区域蔓延。
	if known*2 < len(window) {
		return false
	}
	if hasHomopolymer(window, homopolymerRun) {
		return true
	}
	if hasTandemRepeat(window, 3) {
		return true
	}
	return shannonEntropy(window) < minEntropyBits
}

func maskHomopolymers(seq []byte) {
	run := 1
	for i := 1; i <= len(seq); i++ {
		if i < len(seq) && seq[i] == seq[i-1] && seq[i] != 'N' {
			run++
			continue
		}
		if run >= homopolymerRun {
			for j := i - run; j < i; j++ {
				seq[j] = 'N'
			}
		}
		run = 1
	}
}

func hasHomopolymer(window []byte, minRun int) bool {
	if len(window) < minRun {
		return false
	}
	run := 1
	for i := 1; i < len(window); i++ {
		if window[i] == window[i-1] && window[i] != 'N' {
			run++
			if run >= minRun {
				return true
			}
		} else {
			run = 1
		}
	}
	return false
}

// hasTandemRepeat 在窗口中寻找至少 minRepeats 次完整重复的周期序列。
func hasTandemRepeat(window []byte, minRepeats int) bool {
	for period := 2; period <= maxRepeatPeriod; period++ {
		needed := period * minRepeats
		if needed > len(window) {
			continue
		}
		for start := 0; start+needed <= len(window); start++ {
			if tandemAt(window[start:start+needed], period) {
				return true
			}
		}
	}
	return false
}

func tandemAt(segment []byte, period int) bool {
	if segment[0] == 'N' {
		return false
	}
	// 真正的串联重复：所有位置必须与周期开头对齐。
	// 允许末尾是不完整单元（如长度 7、周期 2 时最后的半个单元），
	// 但不允许开头存在与周期错位的碱基。
	for i := 0; i < len(segment); i++ {
		if segment[i] == 'N' {
			return false
		}
		if i >= period && segment[i] != segment[i%period] {
			return false
		}
	}
	return true
}

// repeatMatch 返回窗口内串联重复的周期与起点（窗口内偏移），找不到时 ok=false。
func repeatMatch(window []byte) (period, segStart int, ok bool) {
	for p := 2; p <= maxRepeatPeriod; p++ {
		needed := p * 3
		if needed > len(window) {
			continue
		}
		for start := 0; start+needed <= len(window); start++ {
			if tandemAt(window[start:start+needed], p) {
				return p, start, true
			}
		}
	}
	return 0, 0, false
}

func shannonEntropy(window []byte) float64 {
	var counts [4]int
	known := 0
	for _, base := range window {
		if idx := baseIndex(base); idx >= 0 {
			counts[idx]++
			known++
		}
	}
	if known == 0 {
		return 0
	}
	entropy := 0.0
	for _, count := range counts {
		if count == 0 {
			continue
		}
		p := float64(count) / float64(known)
		entropy -= p * math.Log2(p)
	}
	return entropy
}
