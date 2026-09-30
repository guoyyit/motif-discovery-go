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
	cfg.MaskLowComplexity = true
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
