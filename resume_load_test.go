package motif

import (
	"context"
	"strings"
	"testing"
)

// TestServiceHighLoadCancelResume 需求 6：增大重排次数与初始化轮数后，
// 在重排阶段中途取消并恢复，结果必须与不中断运行逐位一致。
func TestServiceHighLoadCancelResume(t *testing.T) {
	motifBytes := []byte("GCTAGCTA")
	seqs, _, _ := syntheticInput(t, 321, 12, 50, motifBytes)
	cfg := testConfig(8)
	cfg.Shuffles = 80
	cfg.ExpectMotifs = 3
	cfg.InitRounds = 14
	cfg.MaxEMIterations = 300

	baseline, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}

	svc := newTestService(t)
	ctx := context.Background()
	if err := svc.Start(ctx, "heavy", seqs, cfg); err != nil {
		t.Fatal(err)
	}
	waitForCancelPoint(t, svc, "heavy", "shuffle-37", cfg.Shuffles, cfg.ExpectMotifs)
	if err := svc.Cancel("heavy"); err != nil {
		t.Fatal(err)
	}
	snap, err := svc.Wait(ctx, "heavy")
	if err != nil {
		t.Fatal(err)
	}
	if snap.Status != StatusCanceled || snap.Result != nil {
		t.Fatalf("cancel must leave no partial result: status=%s result=%v",
			snap.Status, snap.Result)
	}
	if err := svc.Resume(ctx, "heavy"); err != nil {
		t.Fatal(err)
	}
	snap, err = svc.Wait(ctx, "heavy")
	if err != nil {
		t.Fatal(err)
	}
	if snap.Status != StatusCompleted {
		t.Fatalf("status=%s err=%s", snap.Status, snap.Error)
	}
	if !resultsEqual(baseline, snap.Result) {
		t.Fatal("high-load resumed result differs from uninterrupted baseline")
	}
}

// TestServiceFailureIsolation 需求 6：单个任务失败（非法输入导致 EM
// 不收敛或配置外的运行时错误）不得影响并发的其他任务。
func TestServiceFailureIsolation(t *testing.T) {
	svc := newTestService(t)
	ctx := context.Background()

	goodSeqs, _, _ := syntheticInput(t, 11, 10, 50, []byte("ACGTATGC"))
	goodCfg := testConfig(8)
	goodCfg.Shuffles = 12

	if err := svc.Start(ctx, "good", goodSeqs, goodCfg); err != nil {
		t.Fatal(err)
	}

	// 坏任务：MaxEMIterations=1 必然在发现阶段以 ErrEMNotConverged 失败。
	badSeqs, _, _ := syntheticInput(t, 22, 10, 50, []byte("TGCATGCA"))
	badCfg := testConfig(8)
	badCfg.Shuffles = 12
	badCfg.MaxEMIterations = 1
	if err := svc.Start(ctx, "bad", badSeqs, badCfg); err != nil {
		t.Fatal(err)
	}

	snapBad, err := svc.Wait(ctx, "bad")
	if err != nil {
		t.Fatal(err)
	}
	if snapBad.Status != StatusFailed {
		t.Fatalf("bad task status=%s, want failed", snapBad.Status)
	}
	// 错误以字符串形式持久化在快照中，需能核对出是不收敛结论。
	if !strings.Contains(snapBad.Error, "EM did not converge") {
		t.Fatalf("bad task error should report non-convergence, got %q", snapBad.Error)
	}

	snapGood, err := svc.Wait(ctx, "good")
	if err != nil {
		t.Fatal(err)
	}
	if snapGood.Status != StatusCompleted {
		t.Fatalf("failure of one task must not affect the other: status=%s err=%s",
			snapGood.Status, snapGood.Error)
	}
	expect, _ := Discover(goodSeqs, goodCfg)
	if !resultsEqual(expect, snapGood.Result) {
		t.Fatal("surviving task result differs from standalone run")
	}
}

// TestServiceConcurrentNoCrossRng 需求 6：多任务并发 + 各自较大的
// 重排次数，结果必须与串行独立运行逐位一致（无跨任务随机状态干扰）。
func TestServiceConcurrentNoCrossRng(t *testing.T) {
	svc := newTestService(t)
	ctx := context.Background()

	specs := []struct {
		id    string
		seed  uint64
		motif []byte
	}{
		{"j0", 101, []byte("ACGTATGC")},
		{"j1", 202, []byte("TGCATGCA")},
		{"j2", 303, []byte("GATCGATC")},
	}
	cfgs := make([]Config, len(specs))
	allSeqs := make([][]Sequence, len(specs))
	for i, sp := range specs {
		seqs, _, _ := syntheticInput(t, sp.seed, 10, 50, sp.motif)
		allSeqs[i] = seqs
		cfgs[i] = testConfig(8)
		cfgs[i].Seed = sp.seed
		cfgs[i].Shuffles = 30
		cfgs[i].InitRounds = 8
		if err := svc.Start(ctx, sp.id, seqs, cfgs[i]); err != nil {
			t.Fatal(err)
		}
	}
	for i, sp := range specs {
		snap, err := svc.Wait(ctx, sp.id)
		if err != nil {
			t.Fatal(err)
		}
		if snap.Status != StatusCompleted {
			t.Fatalf("task %s status=%s err=%s", sp.id, snap.Status, snap.Error)
		}
		expect, err := Discover(allSeqs[i], cfgs[i])
		if err != nil {
			t.Fatal(err)
		}
		if !resultsEqual(expect, snap.Result) {
			t.Fatalf("task %s concurrent result differs from standalone", sp.id)
		}
	}
}
