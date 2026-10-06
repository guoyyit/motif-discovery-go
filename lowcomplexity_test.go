package motif

import (
	"strings"
	"testing"
)

func TestDiscoverFiltersSimpleRepeats(t *testing.T) {
	// 每条序列由随机背景 + 一个长简单重复（ATAT...）片段构成。
	// 不做掩码时该重复会富集为“模体”；开启掩码后不得被报告为显著。
	const width = 8
	r := streamRand(888, globalStream(9001, 0))
	bases := []byte("ACGT")
	seqs := make([]Sequence, 20)
	for i := range seqs {
		seq := make([]byte, 0, 60)
		for len(seq) < 40 {
			seq = append(seq, bases[r.intn(4)])
		}
		repeat := strings.Repeat("AT", 10) // 20bp 的二核苷酸重复
		seq = append(seq, repeat...)
		seqs[i] = Sequence{Data: seq}
	}

	cfg := testConfig(width)
	cfg.Shuffles = 40
	cfg.MaskLowComplexity = Bool(true)
	result, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	for _, m := range result.Motifs {
		c := m.Matrix.Consensus()
		if IsLowComplexity([]byte(c)) {
			t.Fatalf("low-complexity consensus reported as motif: %q (p=%g)", c, m.PValue)
		}
		if m.Significant && IsLowComplexity([]byte(c)) {
			t.Fatalf("simple repeat declared significant: %q", c)
		}
	}
}

// TestExplicitMaskingFalsePreservesOptOut 需求 9：显式关闭掩码是
// 调用方的明确选择，该开关必须继续按原语义生效——重复片段不再被掩码，
// 会进入候选窗口（区别于默认/显式开启路径）。
func TestExplicitMaskingFalsePreservesOptOut(t *testing.T) {
	seq := []byte(strings.Repeat("AT", 16))
	masked := MaskLowComplexity(seq) // 直接调用掩码函数始终掩码
	if string(masked) == string(seq) {
		t.Fatal("MaskLowComplexity should mask a pure AT repeat")
	}

	// 通过配置显式关闭：工作序列应保持原样，窗口数量 = n-L+1（无 N 跳过）。
	seqs := []Sequence{{Name: "s", Data: seq}}
	st, err := newDiscoveryState(seqs, Config{
		MotifLength:       8,
		ExpectMotifs:      1,
		Shuffles:          4,
		MaskLowComplexity: Bool(false),
	})
	if err != nil {
		t.Fatal(err)
	}
	if string(st.working[0]) != string(seq) {
		t.Fatalf("explicit masking=false must leave sequence unmasked: %q", st.working[0])
	}
	if got := len(enumerateWindows(st.working, 8, nil)); got != len(seq)-8+1 {
		t.Fatalf("unmasked repeat should yield %d windows, got %d", len(seq)-8+1, got)
	}
}

func TestDiscoverMasksRepeatsBeforeScan(t *testing.T) {
	seq := []byte("ACGTACGT" + strings.Repeat("CA", 12) + "TTGCAACG")
	masked := MaskLowComplexity(seq)
	// 所有位点不得落在被掩码区域
	for _, b := range masked {
		if b != 'N' && b != 'A' && b != 'C' && b != 'G' && b != 'T' {
			t.Fatalf("unexpected base %q", b)
		}
	}
	windows := enumerateWindows([][]byte{masked}, 8, nil)
	if len(windows) != 1 {
		t.Fatalf("expected only the right clean flank window, got %d windows", len(windows))
	}
	if string(windows[0].data) != "TTGCAACG" {
		t.Fatalf("unexpected surviving window: %q", windows[0].data)
	}
}
