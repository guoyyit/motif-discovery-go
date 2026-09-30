package motif

import (
	"math"
	"testing"
)

func TestBuildAndScorePWM(t *testing.T) {
	windows := [][]byte{
		[]byte("ACGT"),
		[]byte("ACGT"),
		[]byte("ACGA"),
	}
	m, err := BuildPWM(windows, DefaultBackground(), 0.1)
	if err != nil {
		t.Fatalf("BuildPWM error: %v", err)
	}
	if m.Width != 4 {
		t.Fatalf("width=%d, want 4", m.Width)
	}
	if m.Consensus() != "ACGT" {
		t.Fatalf("consensus=%q, want ACGT", m.Consensus())
	}

	best := m.Score([]byte("ACGT"))
	worst := m.MinScore()
	if !(best > worst) {
		t.Fatalf("consensus score %v should exceed min score %v", best, worst)
	}
	if math.Abs(m.NormalizedScore([]byte("ACGT"))-1) > 1e-12 {
		t.Fatal("consensus normalized score should be 1")
	}
	if got := m.Threshold(0); math.Abs(got-m.MinScore()) > 1e-12 {
		t.Fatalf("relative threshold 0 = %v, want min %v", got, m.MinScore())
	}
	if got := m.Threshold(1); math.Abs(got-m.MaxScore()) > 1e-12 {
		t.Fatalf("relative threshold 1 = %v, want max %v", got, m.MaxScore())
	}
}

func TestBuildPWMErrors(t *testing.T) {
	if _, err := BuildPWM(nil, DefaultBackground(), 0.1); err == nil {
		t.Fatal("expected error for no windows")
	}
	if _, err := BuildPWM([][]byte{[]byte("AC"), []byte("ACG")}, DefaultBackground(), 0.1); err == nil {
		t.Fatal("expected error for unequal lengths")
	}
	if _, err := BuildPWM([][]byte{[]byte("ACGN")}, DefaultBackground(), 0.1); err == nil {
		t.Fatal("expected error for unknown base")
	}
}

func TestPWMScorePanics(t *testing.T) {
	m, err := BuildPWM([][]byte{[]byte("ACGT")}, DefaultBackground(), 0.1)
	if err != nil {
		t.Fatal(err)
	}
	func() {
		defer func() {
			if recover() == nil {
				t.Fatal("expected panic on wrong window length")
			}
		}()
		m.Score([]byte("AC"))
	}()
	func() {
		defer func() {
			if recover() == nil {
				t.Fatal("expected panic on N base")
			}
		}()
		m.Score([]byte("ACGN"))
	}()
}
