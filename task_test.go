package motif

import (
	"context"
	"sync"
	"testing"
	"time"
)

func testConfig(perms, rounds int) Config {
	cfg := DefaultConfig()
	cfg.MotifLength = 10
	cfg.InitRounds = rounds
	cfg.Significance.Permutations = perms
	cfg.Filter.Enabled = false
	return cfg
}

func TestReproducibility(t *testing.T) {
	rng := NewRNG(2024)
	seqs := generateSequences(rng, 24, 100)
	rng2 := NewRNG(20240)
	seqs = generateSequences(rng2, 24, 100)
	implant(seqs, "CGTAGCTAGC", NewRNG(20241), 0.85)

	cfg := testConfig(40, 16)
	res1, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	res2, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if !resultsEqual(res1, res2) {
		t.Fatal("identical seed/input produced different results")
	}

	cfg2 := cfg
	cfg2.Seed = cfg.Seed + 1
	res3, err := Discover(seqs, cfg2)
	if err != nil {
		t.Fatal(err)
	}
	if resultsEqual(res1, res3) {
		t.Fatal("different seeds unexpectedly produced identical results")
	}
}

func resultsEqual(a, b *Result) bool {
	if a.Seed != b.Seed || len(a.Motifs) != len(b.Motifs) {
		return false
	}
	for i := range a.Motifs {
		ma, mb := a.Motifs[i], b.Motifs[i]
		if ma.Consensus != mb.Consensus ||
			ma.Score != mb.Score ||
			ma.PValue != mb.PValue ||
			ma.AdjustedP != mb.AdjustedP ||
			ma.Information != mb.Information {
			return false
		}
		for j := range ma.Sites {
			if ma.Sites[j] != mb.Sites[j] {
				return false
			}
		}
	}
	return true
}

func TestConcurrentTasksIsolation(t *testing.T) {
	rng := NewRNG(31337)
	seqs := generateSequences(rng, 20, 100)
	implant(seqs, "AACCGGTTAA", NewRNG(31338), 0.8)

	const n = 4
	var wg sync.WaitGroup
	results := make([]*Result, n)
	errs := make([]error, n)
	for i := 0; i < n; i++ {
		wg.Add(1)
		go func(idx int) {
			defer wg.Done()
			cfg := testConfig(30, 12)
			cfg.Seed = 900 + uint64(idx)
			r, err := Discover(seqs, cfg)
			results[idx] = r
			errs[idx] = err
		}(i)
	}
	wg.Wait()

	// 同一输入、各自固定种子：并跑两次结果必须一致。
	for i := 0; i < n; i++ {
		if errs[i] != nil {
			t.Fatalf("task %d error: %v", i, errs[i])
		}
	}
	// 再次顺序运行，应与并发结果逐位一致。
	for i := 0; i < n; i++ {
		cfg := testConfig(30, 12)
		cfg.Seed = 900 + uint64(i)
		r, err := Discover(seqs, cfg)
		if err != nil {
			t.Fatal(err)
		}
		if !resultsEqual(r, results[i]) {
			t.Fatalf("task %d: concurrent vs sequential results differ", i)
		}
	}
}

func TestTaskManagerLifecycle(t *testing.T) {
	rng := NewRNG(55)
	seqs := generateSequences(rng, 16, 80)
	implant(seqs, "TTAAGGCCAA", NewRNG(56), 0.8)

	mgr := NewTaskManager()
	cfg := testConfig(20, 8)
	id, err := mgr.Submit(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
	defer cancel()
	status, res, err := mgr.Wait(ctx, id)
	if err != nil {
		t.Fatal(err)
	}
	if status != StatusDone || res == nil {
		t.Fatalf("unexpected terminal state: %s res=%v", status, res)
	}
}
