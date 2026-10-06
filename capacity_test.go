package motif

import (
	"context"
	"errors"
	"strings"
	"testing"
	"time"
)

// TestConcurrentHighLoadResume 多任务并发、增大重排次数与初始化轮数后，
// 中断—恢复的结果必须与不中断运行逐位一致，且任务间无随机状态干扰。
func TestConcurrentHighLoadResume(t *testing.T) {
	svc := newTestService(t)
	ctx := context.Background()

	cfg := testConfig(8)
	cfg.ExpectMotifs = 2
	cfg.InitRounds = 12
	cfg.Shuffles = 40

	seqsA, _, _ := syntheticInput(t, 501, 10, 50, []byte("ACGTATGC"))
	seqsB, _, _ := syntheticInput(t, 502, 10, 50, []byte("GCTAGCTA"))
	seqsC, _, _ := syntheticInput(t, 503, 10, 50, []byte("TTGCAACG"))

	// 不中断基线。
	baseA, err := Discover(seqsA, cfg)
	if err != nil {
		t.Fatal(err)
	}
	baseB, err := Discover(seqsB, cfg)
	if err != nil {
		t.Fatal(err)
	}
	baseC, err := Discover(seqsC, cfg)
	if err != nil {
		t.Fatal(err)
	}

	start := time.Now()
	for id, seqs := range map[string][]Sequence{"A": seqsA, "B": seqsB, "C": seqsC} {
		if err := svc.Start(ctx, id, seqs, cfg); err != nil {
			t.Fatal(err)
		}
	}

	// 在重排阶段中段取消 B，A、C 不受影响地跑完。
	waitForCancelPoint(t, svc, "B", "shuffle-15", cfg.Shuffles, cfg.ExpectMotifs)
	if err := svc.Cancel("B"); err != nil {
		t.Fatal(err)
	}
	snapB, err := svc.Wait(ctx, "B")
	if err != nil {
		t.Fatal(err)
	}
	if snapB.Status != StatusCanceled || snapB.Result != nil {
		t.Fatalf("B after cancel: status=%s result=%v", snapB.Status, snapB.Result != nil)
	}

	snapA, err := svc.Wait(ctx, "A")
	if err != nil {
		t.Fatal(err)
	}
	snapC, err := svc.Wait(ctx, "C")
	if err != nil {
		t.Fatal(err)
	}
	if snapA.Status != StatusCompleted || snapC.Status != StatusCompleted {
		t.Fatalf("A=%s C=%s (cancel of B must not affect others)", snapA.Status, snapC.Status)
	}

	// 恢复 B 并完成。
	if err := svc.Resume(ctx, "B"); err != nil {
		t.Fatal(err)
	}
	snapB, err = svc.Wait(ctx, "B")
	if err != nil {
		t.Fatal(err)
	}
	if snapB.Status != StatusCompleted {
		t.Fatalf("B after resume: status=%s err=%s", snapB.Status, snapB.Error)
	}
	t.Logf("3 concurrent tasks (shuffles=%d, initRounds=%d) finished in %v",
		cfg.Shuffles, cfg.InitRounds, time.Since(start))

	if !resultsEqual(baseA, snapA.Result) {
		t.Fatal("A differs from uninterrupted baseline")
	}
	if !resultsEqual(baseB, snapB.Result) {
		t.Fatal("resumed B differs from uninterrupted baseline")
	}
	if !resultsEqual(baseC, snapC.Result) {
		t.Fatal("C differs from uninterrupted baseline")
	}
}

// TestFailingTaskIsolation 单个任务失败或取消不得影响其他并发任务。
func TestFailingTaskIsolation(t *testing.T) {
	svc := newTestService(t)
	ctx := context.Background()

	good := testConfig(8)
	good.Shuffles = 12
	seqsOK, _, _ := syntheticInput(t, 601, 10, 50, []byte("ACGTATGC"))
	baseline, err := Discover(seqsOK, good)
	if err != nil {
		t.Fatal(err)
	}

	// 注定失败的任务：迭代上限为 1，EM 必然不收敛。
	bad := good
	bad.MaxIterations = 1

	if err := svc.Start(ctx, "good", seqsOK, good); err != nil {
		t.Fatal(err)
	}
	if err := svc.Start(ctx, "bad", seqsOK, bad); err != nil {
		t.Fatal(err)
	}

	snapBad, err := svc.Wait(ctx, "bad")
	if err != nil {
		t.Fatal(err)
	}
	if snapBad.Status != StatusFailed {
		t.Fatalf("bad task status=%s, want failed", snapBad.Status)
	}
	if !strings.Contains(snapBad.Error, ErrNotConverged.Error()) &&
		!strings.Contains(snapBad.Error, "converge") {
		t.Fatalf("bad task error should mention non-convergence: %q", snapBad.Error)
	}
	if snapBad.Result != nil {
		t.Fatal("failed task must not expose a partial result")
	}

	snapGood, err := svc.Wait(ctx, "good")
	if err != nil {
		t.Fatal(err)
	}
	if snapGood.Status != StatusCompleted {
		t.Fatalf("good task affected by failing task: status=%s err=%s",
			snapGood.Status, snapGood.Error)
	}
	if !resultsEqual(baseline, snapGood.Result) {
		t.Fatal("good task result differs from standalone baseline")
	}

	// 失败任务的检查点保持 failed 状态，可跨实例读取。
	snap, err := svc.Get("bad")
	if err != nil {
		t.Fatal(err)
	}
	if snap.Status != StatusFailed {
		t.Fatalf("persisted bad task status=%s", snap.Status)
	}
}

// TestCanceledTaskNoPartialState 取消后的任务可恢复出与基线一致的结果，
// 且恢复过程不产生重复计算（零分布条数恰好等于 Shuffles）。
func TestCanceledTaskNoPartialState(t *testing.T) {
	svc := newTestService(t)
	ctx := context.Background()

	cfg := testConfig(8)
	cfg.Shuffles = 30
	seqs, _, _ := syntheticInput(t, 701, 10, 50, []byte("GCTAGCTA"))
	baseline, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}

	if err := svc.Start(ctx, "t", seqs, cfg); err != nil {
		t.Fatal(err)
	}
	waitForCancelPoint(t, svc, "t", "shuffle-10", cfg.Shuffles, cfg.ExpectMotifs)
	if err := svc.Cancel("t"); err != nil {
		t.Fatal(err)
	}
	if _, err := svc.Wait(ctx, "t"); err != nil {
		t.Fatal(err)
	}
	if err := svc.Resume(ctx, "t"); err != nil {
		t.Fatal(err)
	}
	snap, err := svc.Wait(ctx, "t")
	if err != nil {
		t.Fatal(err)
	}
	if snap.Status != StatusCompleted {
		t.Fatalf("status=%s err=%s", snap.Status, snap.Error)
	}
	if !resultsEqual(baseline, snap.Result) {
		t.Fatal("resumed result differs from baseline")
	}

	// 检查点中的零分布恰好累积 Shuffles 条，无重复计算。
	cp, err := loadCheckpoint(svc.dir, "t")
	if err != nil {
		t.Fatal(err)
	}
	if cp.ShuffleDone != cfg.Shuffles {
		t.Fatalf("shuffle_done=%d, want %d (duplicate or missing units)",
			cp.ShuffleDone, cfg.Shuffles)
	}
	for k, nulls := range cp.NullScores {
		if len(nulls) != cfg.Shuffles {
			t.Fatalf("null scores for motif %d: %d entries, want %d",
				k, len(nulls), cfg.Shuffles)
		}
	}
}

// TestDiscoverInvalidConfigError 运行期配置错误必须可通过 errors 链检查。
func TestDiscoverInvalidConfigError(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 702, 6, 40, nil)
	cfg := testConfig(8)
	cfg.MaxIterations = 1
	_, err := Discover(seqs, cfg)
	if !errors.Is(err, ErrNotConverged) {
		t.Fatalf("expected ErrNotConverged, got %v", err)
	}
}
