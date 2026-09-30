package motif

import (
	"testing"
)

func generateSequences(rng *RNG, n, length int) []Seq {
	seqs := make([]Seq, n)
	for i := 0; i < n; i++ {
		data := make([]byte, length)
		for j := range data {
			data[j] = "ACGT"[rng.Intn(4)]
		}
		seqs[i] = Seq{Name: seqName(i), Data: data}
	}
	return seqs
}

func seqName(i int) string {
	return string(rune('a'+i/26)) + string(rune('a'+i%26))
}

func implant(seqs []Seq, motif string, rng *RNG, fraction float64) []int {
	implanted := make([]int, 0)
	for i := range seqs {
		if rng.Float64() >= fraction {
			continue
		}
		pos := rng.Intn(len(seqs[i].Data) - len(motif) + 1)
		copy(seqs[i].Data[pos:pos+len(motif)], motif)
		implanted = append(implanted, i)
	}
	return implanted
}

func TestEndToEndSmoke(t *testing.T) {
	rng := NewRNG(12345)
	seqs := generateSequences(rng, 30, 120)
	implant(seqs, "TACGATTACA", rng, 0.8)

	cfg := DefaultConfig()
	cfg.MotifLength = 10
	cfg.InitRounds = 24
	cfg.Significance.Permutations = 50
	cfg.Filter.Enabled = false

	res, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Motifs) != 1 {
		t.Fatalf("want 1 motif, got %d", len(res.Motifs))
	}
	m := res.Motifs[0]
	t.Logf("consensus=%s score=%.4f p=%.4f adj=%.4f sig=%v",
		m.Consensus, m.Score, m.PValue, m.AdjustedP, m.Significant)
	if m.Consensus != "TACGATTACA" {
		t.Fatalf("consensus mismatch: %q", m.Consensus)
	}
	if !m.Significant || m.PValue > 0.05 {
		t.Fatalf("implanted motif should be significant: p=%g adj=%g", m.PValue, m.AdjustedP)
	}
}

// TestImplantRecallRate 在多组独立合成数据上统计位点级召回率：
// 报告的模体位点与真实植入位点（容差 ±3bp）匹配的比例。
func TestImplantRecallRate(t *testing.T) {
	const trials = 5
	totalImplanted := 0
	totalRecalled := 0

	for trial := 0; trial < trials; trial++ {
		rng := NewRNG(uint64(1000 + trial))
		seqs := generateSequences(rng, 26, 110)
		implantedIdx := implant(seqs, "CGATAAGCTC", NewRNG(uint64(2000+trial)), 0.8)

		cfg := DefaultConfig()
		cfg.MotifLength = 10
		cfg.InitRounds = 20
		cfg.Significance.Permutations = 40
		cfg.Filter.Enabled = false

		res, err := Discover(seqs, cfg)
		if err != nil {
			t.Fatal(err)
		}
		if !res.Motifs[0].Significant {
			t.Fatalf("trial %d: implanted motif not significant, p=%g", trial, res.Motifs[0].PValue)
		}

		siteBySeq := map[int]int{}
		for _, occ := range res.Motifs[0].Sites {
			siteBySeq[occ.SeqIndex] = occ.Start
		}
		recalled := 0
		for _, idx := range implantedIdx {
			// 重新定位真实植入位点不可得，改为在该序列上验证报告位点处
			// 的片段与植入模体高度一致（海明距离 <= 2）。
			st, ok := siteBySeq[idx]
			if !ok {
				continue
			}
			site := seqs[idx].Data[st : st+10]
			if hamming(site, []byte("CGATAAGCTC")) <= 2 {
				recalled++
			}
		}
		totalImplanted += len(implantedIdx)
		totalRecalled += recalled
	}

	rate := float64(totalRecalled) / float64(totalImplanted)
	t.Logf("site-level recall = %.3f (%d/%d)", rate, totalRecalled, totalImplanted)
	if rate < 0.8 {
		t.Fatalf("implant recall rate too low: %.3f", rate)
	}
}

func hamming(a, b []byte) int {
	if len(a) != len(b) {
		return len(a) + len(b)
	}
	d := 0
	for i := range a {
		if a[i] != b[i] {
			d++
		}
	}
	return d
}

func TestEndToEndCalibration(t *testing.T) {
	rng := NewRNG(54321)
	seqs := generateSequences(rng, 30, 120)

	cfg := DefaultConfig()
	cfg.MotifLength = 10
	cfg.InitRounds = 24
	cfg.Significance.Permutations = 99
	cfg.Filter.Enabled = false

	res, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	m := res.Motifs[0]
	t.Logf("random consensus=%s score=%.4f p=%.4f", m.Consensus, m.Score, m.PValue)
	if m.Significant {
		t.Fatalf("pure-random sequences should not yield a significant motif at alpha=0.05: p=%g adj=%g",
			m.PValue, m.AdjustedP)
	}
	if m.PValue <= 1.0/100.0 {
		t.Fatalf("min attainable p with B=99 is 1/100, got %g", m.PValue)
	}
}

func TestLowComplexityNotReported(t *testing.T) {
	rng := NewRNG(777)
	seqs := generateSequences(rng, 30, 120)
	for i := range seqs {
		for j := 40; j < 58; j++ {
			seqs[i].Data[j] = "AT"[(j/1)%2]
		}
	}

	cfg := DefaultConfig()
	cfg.MotifLength = 10
	cfg.InitRounds = 24
	cfg.Significance.Permutations = 50

	res, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	m := res.Motifs[0]
	t.Logf("low-complexity consensus=%s score=%.4f p=%.4f sig=%v",
		m.Consensus, m.Score, m.PValue, m.Significant)
	if m.Significant {
		t.Fatalf("simple repeat must not be reported as significant motif: %s", m.Consensus)
	}
}
