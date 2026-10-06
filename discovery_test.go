package motif

import "testing"

func TestDiscoverImplantRecall(t *testing.T) {
	motif := []byte("ACGTATGC")
	seqs, _, implantedAt := syntheticInput(t, 123, 20, 60, motif)
	cfg := testConfig(8)

	result, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(result.Motifs) != 1 {
		t.Fatalf("got %d motifs, want 1", len(result.Motifs))
	}
	m := result.Motifs[0]
	if m.Matrix.Consensus() != string(motif) {
		t.Fatalf("consensus=%q, want %q", m.Matrix.Consensus(), motif)
	}
	if !m.Significant {
		t.Fatalf("implanted motif not significant: p=%g adj=%g", m.PValue, m.AdjustedPValue)
	}

	// 检查位点召回：每条序列的植入位置必须被报告。
	hit := make([]bool, len(seqs))
	for _, site := range m.Sites {
		if site.Offset == implantedAt[site.SeqIndex] {
			hit[site.SeqIndex] = true
		}
	}
	recovered := 0
	for _, ok := range hit {
		if ok {
			recovered++
		}
	}
	if recovered < 18 {
		t.Fatalf("recovered %d/20 implanted sites", recovered)
	}
}

func TestDiscoverReproducible(t *testing.T) {
	motif := []byte("TGCATGCA")
	seqs, _, _ := syntheticInput(t, 55, 15, 50, motif)
	cfg := testConfig(8)

	r1, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	r2, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if !resultsEqual(r1, r2) {
		t.Fatal("same seed and input must produce identical results")
	}

	// 零值种子同样必须可复现
	cfg2 := cfg
	cfg2.Seed = 0
	r3, err := Discover(seqs, cfg2)
	if err != nil {
		t.Fatal(err)
	}
	r4, err := Discover(seqs, cfg2)
	if err != nil {
		t.Fatal(err)
	}
	if !resultsEqual(r3, r4) {
		t.Fatal("zero seed run is not reproducible")
	}
}

func TestDiscoverValidatesConfig(t *testing.T) {
	seqs := []Sequence{{Name: "s", Data: []byte("ACGTACGTAC")}}
	bad := testConfig(8)
	bad.MotifLength = 2
	if _, err := Discover(seqs, bad); err == nil {
		t.Fatal("expected error for motif length < 4")
	}
	bad = testConfig(8)
	bad.ScoreThreshold = Float64(2)
	if _, err := Discover(seqs, bad); err == nil {
		t.Fatal("expected error for threshold > 1")
	}
}

func resultsEqual(a, b *Result) bool {
	if len(a.Motifs) != len(b.Motifs) {
		return false
	}
	for i := range a.Motifs {
		ma, mb := a.Motifs[i], b.Motifs[i]
		if ma.Matrix.Consensus() != mb.Matrix.Consensus() ||
			ma.PValue != mb.PValue || ma.AdjustedPValue != mb.AdjustedPValue ||
			ma.Significant != mb.Significant || len(ma.Sites) != len(mb.Sites) ||
			ma.Statistic != mb.Statistic || ma.NullCount != mb.NullCount ||
			ma.Information != mb.Information {
			return false
		}
		for j := range ma.Sites {
			if ma.Sites[j] != mb.Sites[j] {
				return false
			}
		}
	}
	return true
}
