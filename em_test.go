package motif

import "testing"

func TestRunEMRecoversImplant(t *testing.T) {
	const width = 8
	const nSeq = 20
	const seqLen = 60

	r := streamRand(123, 0)
	bases := []byte("ACGT")
	motif := []byte("ACGTATGC")

	seqs := make([][]byte, nSeq)
	implanted := make([]int, nSeq)
	for i := 0; i < nSeq; i++ {
		seq := make([]byte, seqLen)
		for j := range seq {
			seq[j] = bases[r.intn(4)]
		}
		pos := r.intn(seqLen - width + 1)
		copy(seq[pos:pos+width], motif)
		seqs[i] = seq
		implanted[i] = pos
	}

	windows := enumerateWindows(seqs, width, nil)
	model, err := fitOneMotif(windows, width, DefaultBackground(), 0.1, 1e-4, 10, 1000,
		999, globalStream(100, 0), "observed")
	if err != nil {
		t.Fatal(err)
	}

	hits := 0
	for i, site := range model.sites {
		if site.SeqIndex == i && site.Offset == implanted[i] {
			hits++
		}
	}
	if hits < 18 {
		t.Fatalf("EM recovered %d/%d implanted sites, want >= 18", hits, nSeq)
	}
	if model.matrix.Consensus() != string(motif) {
		t.Fatalf("consensus=%q, want %q", model.matrix.Consensus(), motif)
	}
}

func TestEnumerateWindowsSkipsN(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACGTACGT"),
		[]byte("ACNTACGT"),
	}
	windows := enumerateWindows(seqs, 4, nil)
	// 序列1：offset 0..2 的窗口都包含 N，只剩 offset 3 (TACG) 与 4 (ACGT)
	count := 0
	for _, w := range windows {
		if w.seqIndex == 1 {
			count++
		}
	}
	if count != 2 {
		t.Fatalf("seq1 windows=%d, want 2: %+v", count, windows)
	}
}
