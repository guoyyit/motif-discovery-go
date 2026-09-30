package motif

import "testing"

func TestEmpiricalPValue(t *testing.T) {
	null := []float64{0.1, 0.2, 0.5, 0.8}
	if got := EmpiricalPValue(0.9, null); got != 1.0/5.0 {
		t.Fatalf("p for above-all = %v", got)
	}
	if got := EmpiricalPValue(0.05, null); got != 5.0/5.0 {
		t.Fatalf("p for below-all = %v", got)
	}
	if got := EmpiricalPValue(0.5, nil); got != 1.0 {
		t.Fatalf("empty null must give p=1, got %v", got)
	}
}

func TestEmpiricalPValueTieCount(t *testing.T) {
	null := []float64{0.1, 0.2, 0.5, 0.8}
	// >= 0.2: 0.2,0.5,0.8 => 3 => p = 4/5.
	if got := EmpiricalPValue(0.2, null); got != 4.0/5.0 {
		t.Fatalf("tie p = %v", got)
	}
}

func TestBHAdjust(t *testing.T) {
	p := []float64{0.01, 0.02, 0.03, 0.5}
	q := BHAdjust(p)
	// BH: sorted rank1 0.01*4/1=0.04, rank2 0.02*4/2=0.04,
	// rank3 0.03*4/3=0.04, rank4 0.5*4/4=0.5 (monotone floor).
	want := []float64{0.04, 0.04, 0.04, 0.5}
	for i := range p {
		if abs(q[i]-want[i]) > 1e-12 {
			t.Fatalf("BH q[%d]=%v want %v (all=%v)", i, q[i], want[i], q)
		}
	}

	// 顺序必须保持。
	q2 := BHAdjust([]float64{0.5, 0.01, 0.02, 0.03})
	if abs(q2[1]-0.04) > 1e-12 || q2[0] != 0.5 {
		t.Fatalf("BH order not preserved: %v", q2)
	}

	if len(BHAdjust(nil)) != 0 {
		t.Fatal("BHAdjust(nil) should be empty")
	}
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func TestMultiMotifOrderedBH(t *testing.T) {
	rng := NewRNG(6060)
	seqs := generateSequences(rng, 40, 160)
	implant(seqs, "AAAAGGTACC", NewRNG(6061), 0.7)
	implant(seqs, "CCTTTTGGAA", NewRNG(6062), 0.7)

	cfg := testConfig(30, 12)
	cfg.ExpectedCount = 2
	res, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Motifs) != 2 {
		t.Fatalf("want 2 motifs, got %d", len(res.Motifs))
	}
	for i, m := range res.Motifs {
		t.Logf("motif %d: %s score=%.3f p=%.4f adj=%.4f sig=%v",
			i, m.Consensus, m.Score, m.PValue, m.AdjustedP, m.Significant)
	}
	// 两个植入模体均应通过校正后阈值。
	sig := 0
	for _, m := range res.Motifs {
		if m.Significant {
			sig++
		}
	}
	if sig < 2 {
		t.Fatalf("expected both implanted motifs significant, got %d", sig)
	}
}
