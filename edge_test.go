package motif

import "testing"

func TestEdgeShortSeqs(t *testing.T) {
	seqs := []Sequence{{Data: []byte("ACG")}}
	cfg := testConfig(8)
	if _, err := Discover(seqs, cfg); err == nil {
		t.Fatal("expected error when sequences are shorter than motif length")
	}
}

func TestEdgeAllMasked(t *testing.T) {
	seqs := []Sequence{
		{Data: []byte("ATATATATATATATATAT")},
		{Data: []byte("ATATATATATATATATAT")},
	}
	cfg := testConfig(8)
	cfg.Shuffles = 5
	cfg.MaskLowComplexity = true
	res, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Motifs) != 0 {
		t.Fatalf("fully masked sequences must yield no motifs, got %d", len(res.Motifs))
	}
}

func TestEdgeEmptyInput(t *testing.T) {
	if _, err := Discover(nil, testConfig(8)); err == nil {
		t.Fatal("expected error for empty input")
	}
}
