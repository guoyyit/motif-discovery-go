package motif

import "math/bits"

// RNG 是由显式种子控制、与 Go map/goroutine 调度无关的确定性随机源。
// 底层为 SplitMix64（状态递增 + 位混合），同一序列的方法调用在任何机器、
// 任何并发调度下都产生完全相同的输出。
type RNG struct {
	state uint64
}

// NewRNG 以单个 64 位种子构造确定性随机源。
func NewRNG(seed uint64) *RNG {
	return &RNG{state: splitMix64(seed + splitMix64Gamma)}
}

// Uint64 返回下一个伪随机 64 位无符号整数。
func (r *RNG) Uint64() uint64 {
	r.state += splitMix64Gamma
	return splitMix64(r.state)
}

// Intn 返回 [0,n) 范围内的伪随机整数。
func (r *RNG) Intn(n int) int {
	if n <= 0 {
		panic("motif: Intn 参数必须为正数")
	}
	// Lemire 有偏拒绝法，保证区间内严格均匀。
	limit := uint64(n)
	hi, lo := bits.Mul64(limit, r.Uint64())
	if lo < limit {
		threshold := -limit % limit
		for lo < threshold {
			hi, lo = bits.Mul64(limit, r.Uint64())
		}
	}
	return int(hi)
}

// Float64 返回 [0,1) 范围内的伪随机浮点数。
func (r *RNG) Float64() float64 {
	return float64(r.Uint64()>>11) / float64(uint64(1)<<53)
}

// Shuffle 使用 Fisher-Yates 算法原地打乱字节序列。
func (r *RNG) Shuffle(b []byte) {
	for i := len(b) - 1; i > 0; i-- {
		j := r.Intn(i + 1)
		b[i], b[j] = b[j], b[i]
	}
}

// Shuffled 返回保留组成但顺序随机化的序列副本。
func (r *RNG) Shuffled(b []byte) []byte {
	out := append([]byte(nil), b...)
	r.Shuffle(out)
	return out
}

// StreamFromSeed 派生子随机源，用于彼此独立的并发任务或随机化轮次。
func StreamFromSeed(seed uint64, stream uint64) *RNG {
	s := splitMix64(seed ^ splitMix64Gamma)
	s ^= splitMix64(stream + splitMix64Gamma)
	return &RNG{state: s}
}

const splitMix64Gamma = 0x9E3779B97F4A7C15

func splitMix64(z uint64) uint64 {
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
	z = (z ^ (z >> 27)) * 0x94D049BB133111EB
	return z ^ (z >> 31)
}
