package motif

import (
	"math/bits"
	mrand "math/rand/v2"
)

// rng 是整个包使用的确定性随机源。
// 基于 math/rand/v2 的 PCG128，所有随机过程（初始化、重排、采样）
// 都必须经由由显式种子派生出的 rng 完成，不依赖全局随机状态。
type rng struct {
	rand   *mrand.Rand
	seed   uint64
	stream uint64
}

// streamRand 由主种子与流编号派生一个独立、确定性的随机源。
// 不同 (seed, stream) 组合产生相互独立的随机序列，保证并发任务互不干扰。
func streamRand(seed uint64, stream uint64) *rng {
	streamed := mixSeed(seed, stream)
	hi, lo := splitSeed(streamed)
	return &rng{rand: mrand.New(mrand.NewPCG(hi, lo)), seed: seed, stream: stream}
}

// 流编号采用两级层级：65536 个高级流（每类用途一个），每个高级流含 2^48 个子流。
// 这样嵌套派生（如 shuffle b 内的 motif k 的 EM round r）也不会冲突。
const substreamSpace = uint64(1) << 48

// globalStream 由高级流编号与子流编号组合出全局唯一流编号。
func globalStream(high, sub uint64) uint64 {
	return high*substreamSpace + sub
}

// mixSeed 将流编号以确定性的无冲突方式混入主种子。
func mixSeed(seed, stream uint64) uint64 {
	const c = 0x9E3779B97F4A7C15
	x := seed ^ (stream*c + 0xD1B54A32D192ED03)
	x ^= x >> 30
	x *= 0xBF58476D1CE4E5B9
	x ^= x >> 27
	x *= 0x94D049BB133111EB
	x ^= x >> 31
	return x
}

// splitSeed 将一个 64 位种子展开为 PCG 需要的两个 64 位状态字。
func splitSeed(seed uint64) (uint64, uint64) {
	if seed == 0 {
		// 零状态会让 PCG 停在固定点，这里换成非零常数。
		seed = 0x9E3779B97F4A7C15
	}
	hi := bits.RotateLeft64(seed^0xD1B54A32D192ED03, 17)
	lo := bits.ReverseBytes64(seed) ^ 0xA24BAED4963EE407
	if hi == 0 && lo == 0 {
		lo = 0x9E3779B97F4A7C15
	}
	return hi, lo
}

func (r *rng) uint64() uint64 { return r.rand.Uint64() }

// intn 返回 [0,n) 内的确定性随机整数。
func (r *rng) intn(n int) int {
	if n <= 0 {
		panic("rng: non-positive bound")
	}
	return int(r.rand.Uint64N(uint64(n)))
}

// float64 返回 [0,1) 内的确定性随机浮点数。
func (r *rng) float64() float64 { return r.rand.Float64() }

// perm 返回 [0,n) 的一个确定性随机排列。
func (r *rng) perm(n int) []int { return r.rand.Perm(n) }
