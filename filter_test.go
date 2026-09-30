package motif

import "testing"

func TestFilterMasksMonomerRun(t *testing.T) {
	f := NewComplexityFilter(DefaultFilterConfig())
	seq := []byte("ACGTACGTACGTAAAAAAAAAAAAACGTACGT")
	regions := f.MaskRegions(seq)
	if len(regions) == 0 {
		t.Fatal("expected masked regions for A run")
	}
	foundRun := false
	for _, r := range regions {
		if r.Kind == "monomer-run" && r.Start == 12 && r.End >= 24 {
			foundRun = true
		}
	}
	if !foundRun {
		t.Fatalf("monomer run not detected: %+v", regions)
	}
	if !f.IsMasked(regions, 15) {
		t.Fatal("position 15 should be masked")
	}
}

func TestFilterMasksDinucleotideRepeat(t *testing.T) {
	f := NewComplexityFilter(DefaultFilterConfig())
	seq := []byte("GGC" + repeatUnit("AT", 8) + "GGC")
	regions := f.MaskRegions(seq)
	found := false
	for _, r := range regions {
		if r.Kind == "tandem-repeat" {
			found = true
		}
	}
	if !found {
		t.Fatalf("dinucleotide repeat not detected: %+v", regions)
	}
}

func repeatUnit(unit string, times int) string {
	out := make([]byte, 0, len(unit)*times)
	for i := 0; i < times; i++ {
		out = append(out, unit...)
	}
	return string(out)
}

func TestFilterAllowedStarts(t *testing.T) {
	f := NewComplexityFilter(DefaultFilterConfig())
	seq := make([]byte, 60)
	for i := range seq {
		seq[i] = "ACGT"[i%4]
	}
	allowed := f.AllowedStarts(seq, 8)
	if len(allowed) != 0 {
		t.Fatalf("periodic sequence should have no allowed starts, got %d", len(allowed))
	}

	cfg := DefaultFilterConfig()
	cfg.Enabled = false
	fOff := NewComplexityFilter(cfg)
	if len(fOff.AllowedStarts(seq, 8)) != 60-8+1 {
		t.Fatal("disabled filter should allow all starts")
	}
}

func TestFilterSitePasses(t *testing.T) {
	f := NewComplexityFilter(DefaultFilterConfig())
	blocked := [][]byte{
		[]byte("AAAAAAAAAA"),
		[]byte("ATATATATAT"),
		[]byte("ACACACACAC"),
		[]byte("ACGTACGTAC"),
		[]byte("ACGACGACGA"),
	}
	for _, s := range blocked {
		if f.SitePasses(s) {
			t.Errorf("site %q should be blocked", s)
		}
	}
	good := [][]byte{
		[]byte("ACGTGGTACA"),
		[]byte("TAGCCGATGC"),
	}
	for _, s := range good {
		if !f.SitePasses(s) {
			t.Errorf("site %q should pass", s)
		}
	}
}
