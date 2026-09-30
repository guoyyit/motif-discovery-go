package motif

import "testing"

func TestEstimateBackground(t *testing.T) {
	seqs := []Sequence{
		{Data: []byte("AACN")},
		{Data: []byte("gt")},
	}
	bg := EstimateBackground(seqs)
	total := 0.0
	for _, p := range bg.Freq {
		if p < 0 {
			t.Fatalf("negative frequency: %v", bg.Freq)
		}
		total += p
	}
	if total < 1-1e-12 || total > 1+1e-12 {
		t.Fatalf("frequencies sum to %v, want 1", total)
	}
	if bg.Prob('N') != 0 {
		t.Fatalf("N probability should be 0, got %v", bg.Prob('N'))
	}
	if err := bg.Validate(); err != nil {
		t.Fatalf("valid background rejected: %v", err)
	}
}

func TestDefaultBackground(t *testing.T) {
	bg := DefaultBackground()
	for i, p := range bg.Freq {
		if p != 0.25 {
			t.Fatalf("freq[%d]=%v, want 0.25", i, p)
		}
	}
	if err := bg.Validate(); err != nil {
		t.Fatalf("uniform background rejected: %v", err)
	}
}

func TestBackgroundValidateError(t *testing.T) {
	bad := Background{Freq: [4]float64{0.5, 0.5, 0.5, 0.5}}
	if err := bad.Validate(); err == nil {
		t.Fatal("expected validation error for frequencies summing to 2")
	}
}
