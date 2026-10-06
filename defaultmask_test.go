package motif

import (
	"math"
	"strings"
	"testing"
)

// atRepeatInput 生成随机背景中嵌入 AT 串联重复的序列集合。
func atRepeatInput(t *testing.T, seed uint64, nSeq, flank int) []Sequence {
	t.Helper()
	r := streamRand(seed, globalStream(9100, 0))
	bases := []byte("ACGT")
	seqs := make([]Sequence, nSeq)
	for i := range seqs {
		seq := make([]byte, 0, 2*flank+24)
		for len(seq) < flank {
			seq = append(seq, bases[r.intn(4)])
		}
		seq = append(seq, []byte(strings.Repeat("AT", 12))...) // 24bp 二核苷酸重复
		for len(seq) < 2*flank+24 {
			seq = append(seq, bases[r.intn(4)])
		}
		seqs[i] = Sequence{Name: string(rune('a' + i)), Data: seq}
	}
	return seqs
}

// checkNoSignificantRepeat 验证简单重复不被报告为显著模体，
// 且每个模体都保留可核对的统计量。
func checkNoSignificantRepeat(t *testing.T, result *Result) {
	t.Helper()
	for _, m := range result.Motifs {
		consensus := m.Matrix.Consensus()
		if IsLowComplexity([]byte(consensus)) {
			t.Fatalf("low-complexity consensus reported as motif: %q", consensus)
		}
		if m.Significant {
			t.Fatalf("motif %q declared significant on simple-repeat input (p=%g adj=%g)",
				consensus, m.PValue, m.AdjustedPValue)
		}
		// 未达显著时必须保留可核对的统计量。
		if m.PValue <= 0 || m.PValue > 1 || math.IsNaN(m.PValue) {
			t.Fatalf("p-value not checkable: %g", m.PValue)
		}
		if m.AdjustedPValue < m.PValue-1e-12 || m.AdjustedPValue > 1 || math.IsNaN(m.AdjustedPValue) {
			t.Fatalf("adjusted p-value not checkable: %g (raw %g)", m.AdjustedPValue, m.PValue)
		}
		if math.IsNaN(m.Information) || math.IsInf(m.Information, 0) {
			t.Fatalf("information content not checkable: %g", m.Information)
		}
	}
}

// TestDefaultConfigRejectsSimpleRepeats 默认配置（不设置任何掩码开关）
// 不得把 AT 串联重复报告为显著模体——防护是默认行为而非可选开关。
func TestDefaultConfigRejectsSimpleRepeats(t *testing.T) {
	seqs := atRepeatInput(t, 2026, 8, 20)
	result, err := Discover(seqs, Config{Seed: 42, Shuffles: 60})
	if err != nil {
		t.Fatal(err)
	}
	checkNoSignificantRepeat(t, result)
}

// TestDefaultMatchesExplicitMask 默认路径与显式开启防护（MaskLowComplexity=true）
// 的结论必须逐位一致。
func TestDefaultMatchesExplicitMask(t *testing.T) {
	seqs := atRepeatInput(t, 2027, 8, 20)

	def := Config{Seed: 7, MotifLength: 8, ExpectMotifs: 2, InitRounds: 5, Shuffles: 25}
	explicit := def
	explicit.MaskLowComplexity = true

	rDefault, err := Discover(seqs, def)
	if err != nil {
		t.Fatal(err)
	}
	rExplicit, err := Discover(seqs, explicit)
	if err != nil {
		t.Fatal(err)
	}
	if !resultsEqual(rDefault, rExplicit) {
		t.Fatal("default config and explicit MaskLowComplexity=true must agree bitwise")
	}
	checkNoSignificantRepeat(t, rDefault)
}

// TestDefaultConfigPureRepeatInput 纯 AT 重复输入在默认配置下
// 被完全掩码，不得产生任何模体，也不得报错。
func TestDefaultConfigPureRepeatInput(t *testing.T) {
	seqs := []Sequence{
		{Name: "r1", Data: []byte(strings.Repeat("AT", 30))},
		{Name: "r2", Data: []byte(strings.Repeat("AT", 30))},
	}
	result, err := Discover(seqs, Config{Seed: 1, MotifLength: 8, Shuffles: 5})
	if err != nil {
		t.Fatal(err)
	}
	if len(result.Motifs) != 0 {
		t.Fatalf("pure repeat input must yield no motifs by default, got %d", len(result.Motifs))
	}
}
