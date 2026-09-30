package motif

import (
	"context"
	"path/filepath"
	"testing"
	"time"
)

func TestPauseResumeMatchesUninterrupted(t *testing.T) {
	rng := NewRNG(8080)
	seqs := generateSequences(rng, 24, 100)
	implant(seqs, "ATCGATCGAA", NewRNG(8081), 0.85)

	cfg := testConfig(60, 20)

	// 不中断的参考结果。
	reference, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}

	// 通过 TaskManager 提交，立即请求暂停，等待暂停落盘，再恢复。
	mgr := NewTaskManager()
	id, err := mgr.Submit(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := mgr.Pause(id); err != nil {
		t.Fatal(err)
	}

	deadline := time.Now().Add(10 * time.Second)
	var status TaskStatus
	for time.Now().Before(deadline) {
		status, _, _ = mgr.Get(id)
		if status == StatusPaused {
			break
		}
		time.Sleep(time.Millisecond)
	}
	if status != StatusPaused {
		t.Fatalf("task did not pause in time, status=%s", status)
	}

	// 检查点必须记录了部分进度（存在已完成轮次或重排）。
	cp, err := mgr.Checkpoint(id)
	if err != nil {
		t.Fatal(err)
	}
	cp.mu.Lock()
	progress := len(cp.DoneRounds) + len(cp.NullSamples)
	t.Logf("paused status=%d rounds=%d null=%d motifs=%d", cp.Status, len(cp.DoneRounds), len(cp.NullSamples), len(cp.Motifs))
	cp.mu.Unlock()
	_ = progress

	// 恢复前不得暴露部分最终结果。
	if _, res, _ := mgr.Get(id); res != nil {
		t.Fatal("partial result must not be exposed while paused")
	}

	progressAtPause := func() int {
		cp.mu.Lock()
		defer cp.mu.Unlock()
		return len(cp.DoneRounds) + len(cp.NullSamples)
	}
	if progressAtPause() == 0 {
		if err := mgr.Resume(id); err != nil {
			t.Fatal(err)
		}
		time.Sleep(20 * time.Millisecond)
		if err := mgr.Pause(id); err != nil {
			t.Fatal(err)
		}
		waitFor := func(want TaskStatus) {
			dl := time.Now().Add(10 * time.Second)
			for time.Now().Before(dl) {
				if s, _, _ := mgr.Get(id); s == want {
					return
				}
				time.Sleep(time.Millisecond)
			}
		}
		waitFor(StatusPaused)
		if progressAtPause() == 0 {
			t.Fatal("expected checkpoint progress after a brief run")
		}
	}

	if err := mgr.Resume(id); err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
	defer cancel()
	status, resumed, err := mgr.Wait(ctx, id)
	if err != nil {
		t.Fatal(err)
	}
	if status != StatusDone || resumed == nil {
		t.Fatalf("resumed task did not finish: %s", status)
	}
	if !resultsEqual(reference, resumed) {
		t.Fatal("resumed result differs from uninterrupted run")
	}
}

func TestCheckpointSaveLoadRoundTrip(t *testing.T) {
	rng := NewRNG(9090)
	seqs := generateSequences(rng, 20, 90)
	implant(seqs, "GGCCAATTGG", NewRNG(9091), 0.8)

	cfg := testConfig(40, 12)
	cp := newCheckpoint(cfg)
	dir := t.TempDir()

	// 直接跑一小步：提取模体阶段，检查点记录轮次。
	bg := estimateBackground(seqs)
	filter := NewComplexityFilter(cfg.Filter)
	allowed := make([][]int, len(seqs))
	for i, s := range seqs {
		allowed[i] = filter.AllowedStarts(s.Data, cfg.MotifLength)
	}
	_, err := bestOfRounds(seqs, cfg, bg, allowed, makeBlocked(len(seqs)), cp, 0)
	if err != nil {
		t.Fatal(err)
	}
	path, err := cp.Save(dir)
	if err != nil {
		t.Fatal(err)
	}
	if filepath.Base(path) != cp.TaskID+".json" {
		t.Fatalf("unexpected checkpoint path %s", path)
	}

	loaded, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatal(err)
	}
	if loaded.TaskID != cp.TaskID || loaded.Seed != cp.Seed {
		t.Fatal("checkpoint identity mismatch after load")
	}
	if len(loaded.DoneRounds) != cfg.InitRounds {
		t.Fatalf("expected %d done rounds, got %d", cfg.InitRounds, len(loaded.DoneRounds))
	}
	if err := loaded.match(cfg); err != nil {
		t.Fatal(err)
	}

	// 用加载的检查点继续完整流程，结果与不中断一致。
	reference, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	resumed, err := DiscoverWithProgress(seqs, cfg, loaded)
	if err != nil {
		t.Fatal(err)
	}
	if !resultsEqual(reference, resumed) {
		t.Fatal("loaded-checkpoint resume differs from uninterrupted run")
	}
}

func makeBlocked(n int) []map[int]bool {
	b := make([]map[int]bool, n)
	for i := range b {
		b[i] = map[int]bool{}
	}
	return b
}

func TestDiskBackedPauseResume(t *testing.T) {
	rng := NewRNG(4242)
	seqs := generateSequences(rng, 22, 95)
	implant(seqs, "ACGTACGTAA", NewRNG(4243), 0.85)

	cfg := testConfig(50, 16)
	reference, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}

	dir := t.TempDir()
	cp := newCheckpoint(cfg)
	if _, err := cp.Save(dir); err != nil {
		t.Fatal(err)
	}

	mgr1 := NewTaskManager()
	id, err := mgr1.SubmitWithCheckpoint(seqs, cfg, cp)
	if err != nil {
		t.Fatal(err)
	}
	if err := mgr1.Pause(id); err != nil {
		t.Fatal(err)
	}

	deadline := time.Now().Add(10 * time.Second)
	for time.Now().Before(deadline) {
		if s, _, _ := mgr1.Get(id); s == StatusPaused {
			break
		}
		time.Sleep(time.Millisecond)
	}

	// 模拟新进程：从磁盘加载检查点，用新的管理器恢复。
	path := filepath.Join(dir, cp.TaskID+".json")
	loaded, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatal(err)
	}
	mgr2 := NewTaskManager()
	id2, err := mgr2.SubmitResumable(seqs, cfg, loaded)
	if err != nil {
		t.Fatal(err)
	}
	if id2 != id {
		t.Fatalf("task id changed across resume: %s != %s", id2, id)
	}
	ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
	defer cancel()
	status, final, err := mgr2.Wait(ctx, id2)
	if err != nil {
		t.Fatal(err)
	}
	if status != StatusDone {
		t.Fatalf("expected done, got %s", status)
	}
	if !resultsEqual(reference, final) {
		t.Fatal("disk-backed resumed result differs")
	}
}
