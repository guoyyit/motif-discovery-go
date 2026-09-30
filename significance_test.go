package motif

import (
	"math"
	"testing"
)

func TestEmpiricalPValue(t *testing.T) {
	null := []float64{0.1, 0.2, 0.3, 0.4}
	if p := empiricalPValue(1.0, null); p != 1.0/5 {
		t.Fatalf("extreme p=%v, want %v", p, 1.0/5)
	}
	if p := empiricalPValue(0.0, null); p != 1.0 {
		t.Fatalf("weak p=%v, want 1.0", p)
	}
	if p := empiricalPValue(0.3, null); p != 3.0/5 {
		t.Fatalf("middle p=%v, want %v", p, 3.0/5)
	}
	if p := empiricalPValue(0.5, nil); p != 1.0 {
		t.Fatalf("empty null p=%v, want 1.0", p)
	}
}

func TestBenjaminiHochberg(t *testing.T) {
	// 经典示例：p = [0.01, 0.01, 0.5]，m=3
	q := benjaminiHochberg([]float64{0.01, 0.01, 0.5})
	if math.Abs(q[0]-0.015) > 1e-12 || math.Abs(q[1]-0.015) > 1e-12 {
		t.Fatalf("q[0]=%v q[1]=%v, want 0.015", q[0], q[1])
	}
	if math.Abs(q[2]-0.5) > 1e-12 {
		t.Fatalf("q[2]=%v, want 0.5", q[2])
	}
}

func TestBenjaminiHochbergMonotonic(t *testing.T) {
	pv := []float64{0.001, 0.9, 0.005, 0.02, 0.4, 0.01}
	q := benjaminiHochberg(pv)
	for _, v := range q {
		if v < 0 || v > 1 {
			t.Fatalf("adjusted p out of range: %v", q)
		}
	}
	// 按原始 p 值升序检查校正后 q 值必须单调不减
	ranks := rankOrder(pv)
	prev := 0.0
	for _, idx := range ranks {
		if q[idx] < prev-1e-15 {
			t.Fatalf("non-monotonic adjusted p: %v", q)
		}
		prev = q[idx]
	}
}

func rankOrder(pv []float64) []int {
	idx := make([]int, len(pv))
	for i := range idx {
		idx[i] = i
	}
	for i := 1; i < len(idx); i++ {
		for j := i; j > 0 && pv[idx[j-1]] > pv[idx[j]]; j-- {
			idx[j-1], idx[j] = idx[j], idx[j-1]
		}
	}
	return idx
}
