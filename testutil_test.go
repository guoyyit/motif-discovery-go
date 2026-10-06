package motif

import (
	"fmt"
	"testing"
)

// syntheticInput 生成 nSeq 条长度 seqLen 的随机序列；
// implant 非空时在每条序列的随机位置植入该模体。
func syntheticInput(t *testing.T, seed uint64, nSeq, seqLen int, implant []byte) ([]Sequence, [][]byte, []int) {
	t.Helper()
	r := streamRand(seed, globalStream(9000, 0))
	bases := []byte("ACGT")
	seqs := make([]Sequence, nSeq)
	raw := make([][]byte, nSeq)
	implantedAt := make([]int, nSeq)
	width := len(implant)
	for i := 0; i < nSeq; i++ {
		seq := make([]byte, seqLen)
		for j := range seq {
			seq[j] = bases[r.intn(4)]
		}
		if width > 0 {
			pos := r.intn(seqLen - width + 1)
			copy(seq[pos:pos+width], implant)
			implantedAt[i] = pos
		}
		raw[i] = seq
		seqs[i] = Sequence{Name: fmt.Sprintf("seq%d", i), Data: seq}
	}
	return seqs, raw, implantedAt
}

func testConfig(width int) Config {
	return Config{
		MotifLength:       width,
		ExpectMotifs:      1,
		InitRounds:        10,
		ConvergenceTol:    1e-4,
		ScoreThreshold:    Float64(0.8),
		SignificanceAlpha: Float64(0.05),
		Shuffles:          50,
		Pseudocount:       0.5,
	}
}
