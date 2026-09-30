package motif

import "testing"

func TestStreamRandDeterminism(t *testing.T) {
	first := streamRand(42, 7)
	second := streamRand(42, 7)
	other := streamRand(42, 8)

	same := true
	different := false
	for i := 0; i < 100; i++ {
		a := first.uint64()
		b := second.uint64()
		c := other.uint64()
		if a != b {
			same = false
		}
		if a != c {
			different = true
		}
	}
	if !same {
		t.Fatal("same (seed, stream) must produce identical sequences")
	}
	if !different {
		t.Fatal("different streams should produce different sequences")
	}
}

func TestStreamRandPerm(t *testing.T) {
	r1 := streamRand(1, 0)
	r2 := streamRand(1, 0)
	p1 := r1.perm(20)
	p2 := r2.perm(20)
	if len(p1) != 20 {
		t.Fatalf("perm length %d, want 20", len(p1))
	}
	seen := make(map[int]bool)
	for _, v := range p1 {
		if v < 0 || v >= 20 || seen[v] {
			t.Fatalf("invalid permutation: %v", p1)
		}
		seen[v] = true
	}
	for i := range p1 {
		if p1[i] != p2[i] {
			t.Fatal("permutation not reproducible")
		}
	}
}

func TestStreamRandZeroSeed(t *testing.T) {
	r := streamRand(0, 0)
	_ = r.uint64()
}
