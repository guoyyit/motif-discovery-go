package motif

import "testing"

func TestRNGDeterministic(t *testing.T) {
	a := NewRNG(42)
	b := NewRNG(42)
	for i := 0; i < 1000; i++ {
		if a.Uint64() != b.Uint64() {
			t.Fatalf("Uint64 diverged at %d", i)
		}
	}
	c := NewRNG(43)
	if c.Uint64() == NewRNG(42).Uint64() {
		t.Fatal("different seeds produced same first value")
	}
}

func TestRNGIntnRange(t *testing.T) {
	r := NewRNG(7)
	counts := make([]int, 4)
	for i := 0; i < 4000; i++ {
		v := r.Intn(4)
		if v < 0 || v >= 4 {
			t.Fatalf("out of range: %d", v)
		}
		counts[v]++
	}
	for b, n := range counts {
		if n < 800 || n > 1200 {
			t.Fatalf("bucket %d suspicious count %d", b, n)
		}
	}
}

func TestRNGShuffleComposition(t *testing.T) {
	orig := []byte("AAAACCCGGT")
	r1 := NewRNG(99)
	out := r1.Shuffled(orig)
	if string(orig) != "AAAACCCGGT" {
		t.Fatalf("Shuffled mutated input: %q", orig)
	}
	if !sameComposition(orig, out) {
		t.Fatalf("composition changed: %q -> %q", orig, out)
	}
	a := NewRNG(99).Shuffled(orig)
	b := NewRNG(99).Shuffled(orig)
	if string(a) != string(b) {
		t.Fatal("shuffle not reproducible")
	}
	if StreamFromSeed(1, 0).Uint64() == StreamFromSeed(1, 1).Uint64() {
		t.Fatal("streams should differ")
	}
}

func sameComposition(a, b []byte) bool {
	if len(a) != len(b) {
		return false
	}
	counts := [256]int{}
	for _, c := range a {
		counts[c]++
	}
	for _, c := range b {
		counts[c]--
		if counts[c] < 0 {
			return false
		}
	}
	return true
}
