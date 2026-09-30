package motif

import "testing"

func TestShuffleBases(t *testing.T) {
	original := []byte("ACGTACGTNNACGT")

	out1 := append([]byte(nil), original...)
	shuffleBases(out1, streamRand(7, 0))
	out2 := append([]byte(nil), original...)
	shuffleBases(out2, streamRand(7, 0))

	if string(out1) != string(out2) {
		t.Fatal("same seed must give identical shuffles")
	}

	// N 必须留在原位
	for i, base := range original {
		wasN := base == 'N'
		nowN := out1[i] == 'N'
		if wasN != nowN {
			t.Fatalf("N positions changed at %d: %q -> %q", i, original, out1)
		}
	}

	// 碱基组成必须保持
	if composition(original) != composition(out1) {
		t.Fatalf("composition changed: %q -> %q", original, out1)
	}
}

func composition(seq []byte) [4]int {
	var counts [4]int
	for _, base := range seq {
		if idx := baseIndex(base); idx >= 0 {
			counts[idx]++
		}
	}
	return counts
}

func TestShuffleAllPreservesInput(t *testing.T) {
	seqs := [][]byte{[]byte("AAAA"), []byte("CCCC")}
	out := shuffleAll(seqs, streamRand(1, 0))
	out[0][0] = 'T'
	if string(seqs[0]) != "AAAA" {
		t.Fatal("shuffleAll must not alias input slices")
	}
}
