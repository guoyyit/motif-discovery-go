package motif

import (
	"math"
	"testing"
)

func TestBuildPWMAndScore(t *testing.T) {
	sites := [][]byte{
		[]byte("ACGT"),
		[]byte("ACGT"),
		[]byte("ACGT"),
		[]byte("ACGA"),
	}
	bg := UniformBackground()
	m, err := BuildPWM(sites, bg, DefaultPWMConfig())
	if err != nil {
		t.Fatal(err)
	}
	if m.Length != 4 || m.Consensus != "ACGT" {
		t.Fatalf("length/consensus mismatch: %d %q", m.Length, m.Consensus)
	}
	best := m.Score([]byte("ACGT"))
	worst := m.MinScore()
	if !(best > worst) {
		t.Fatalf("consensus score %v should exceed worst %v", best, worst)
	}
	ns := m.NormalizedScore([]byte("ACGT"))
	if math.Abs(ns-1) > 1e-12 {
		t.Fatalf("normalized consensus score = %v, want 1", ns)
	}
	if !(m.MaxScore() >= best) {
		t.Fatalf("MaxScore %v should be >= consensus score %v", m.MaxScore(), best)
	}
}

func TestPWMScan(t *testing.T) {
	sites := [][]byte{[]byte("AAAA"), []byte("AAAA"), []byte("AAAT")}
	m, err := BuildPWM(sites, UniformBackground(), DefaultPWMConfig())
	if err != nil {
		t.Fatal(err)
	}
	hits := m.Scan([]byte("CCAAAAGG"), 0.99)
	if len(hits) != 1 || hits[0].Start != 2 {
		t.Fatalf("unexpected hits: %+v", hits)
	}
	seqs := []Seq{{Name: "x", Data: []byte("TTAAAATT")}}
	all := m.ScanAll(seqs, 0.99)
	if len(all) != 1 || all[0].SeqIndex != 0 || all[0].Start != 2 {
		t.Fatalf("unexpected scan-all hits: %+v", all)
	}
}

func TestBackgroundValidate(t *testing.T) {
	if err := UniformBackground().Validate(); err != nil {
		t.Fatal(err)
	}
	if err := (Background{0.5, 0.5, 0, 0}).Validate(); err != nil {
		t.Fatal(err)
	}
	if err := (Background{0.5, 0.2, 0.2, 0.2}).Validate(); err == nil {
		t.Fatal("expected invalid background error")
	}
}

func TestBuildPWMErrors(t *testing.T) {
	if _, err := BuildPWM(nil, UniformBackground(), DefaultPWMConfig()); err == nil {
		t.Fatal("expected error on empty sites")
	}
	if _, err := BuildPWM([][]byte{[]byte("AC"), []byte("A")}, UniformBackground(), DefaultPWMConfig()); err == nil {
		t.Fatal("expected error on mismatched lengths")
	}
}
