package motif

import "testing"

func TestIsLowComplexity(t *testing.T) {
	cases := []struct {
		name string
		seq  string
		want bool
	}{
		{"homopolymer", "AAAAAAGT", true},
		{"dinucleotide repeat", "ATATATCG", true},
		{"trinucleotide repeat", "ATGATGATG", true},
		{"low entropy", "AAAAAATT", true},
		{"complex", "ACGTACGT", false},
		{"random-like", "ACGGTTCAG", false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := IsLowComplexity([]byte(tc.seq)); got != tc.want {
				t.Fatalf("IsLowComplexity(%q)=%v, want %v", tc.seq, got, tc.want)
			}
		})
	}
}

func TestMaskLowComplexity(t *testing.T) {
	seq := []byte("GGGGGGGGGACGTACGT")
	masked := MaskLowComplexity(seq)
	if string(masked)[:9] != "NNNNNNNNN" {
		t.Fatalf("repeat run not masked: %q", masked)
	}
	if string(masked)[9:] != "ACGTACGT" {
		t.Fatalf("complex tail wrongly masked: %q", masked)
	}
	// 原始切片不被修改
	if string(seq) != "GGGGGGGGGACGTACGT" {
		t.Fatalf("input slice was modified: %q", seq)
	}
}
