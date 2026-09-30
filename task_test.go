package motif

import (
	"context"
	"path/filepath"
	"testing"
	"time"
)

func spinWait() {
	time.Sleep(time.Millisecond)
}

func newTestService(t *testing.T) *Service {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "state")
	svc, err := NewService(dir)
	if err != nil {
		t.Fatal(err)
	}
	return svc
}

func TestServiceRunToCompletion(t *testing.T) {
	svc := newTestService(t)
	motifBytes := []byte("ACGTATGC")
	seqs, _, _ := syntheticInput(t, 321, 12, 50, motifBytes)
	cfg := testConfig(8)
	cfg.Shuffles = 20

	ctx := context.Background()
	if err := svc.Start(ctx, "job1", seqs, cfg); err != nil {
		t.Fatal(err)
	}
	snap, err := svc.Wait(ctx, "job1")
	if err != nil {
		t.Fatal(err)
	}
	if snap.Status != StatusCompleted {
		t.Fatalf("status=%s err=%s", snap.Status, snap.Error)
	}
	if snap.Result == nil || len(snap.Result.Motifs) != 1 {
		t.Fatalf("unexpected result: %+v", snap.Result)
	}
	if snap.Result.Motifs[0].Matrix.Consensus() != string(motifBytes) {
		t.Fatalf("consensus=%q", snap.Result.Motifs[0].Matrix.Consensus())
	}

	// 完成后再 Start 同 ID 应报错，要求使用 Resume（而 Resume 也应拒绝重复完成）
	if err := svc.Start(ctx, "job1", seqs, cfg); err == nil {
		t.Fatal("duplicate Start should fail")
	}
}

func TestServiceCancelAndResume(t *testing.T) {
	motifBytes := []byte("GCTAGCTA")
	seqs, _, _ := syntheticInput(t, 321, 12, 50, motifBytes)
	cfg := testConfig(8)
	cfg.Shuffles = 24
	cfg.ExpectMotifs = 3

	// 不中断基线
	baseline, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}

	for _, cancelAfter := range []string{"discover-1", "discover", "shuffle-1", "shuffle-12"} {
		name := cancelAfter
		t.Run(name, func(t *testing.T) {
			svc := newTestService(t)
			ctx := context.Background()
			if err := svc.Start(ctx, "t", seqs, cfg); err != nil {
				t.Fatal(err)
			}

			waitForCancelPoint(t, svc, "t", cancelAfter, cfg.Shuffles, cfg.ExpectMotifs)
			if err := svc.Cancel("t"); err != nil {
				t.Fatal(err)
			}
			snap, err := svc.Wait(ctx, "t")
			if err != nil {
				t.Fatal(err)
			}
			if snap.Status != StatusCanceled {
				t.Fatalf("status after cancel=%s err=%s", snap.Status, snap.Error)
			}
			if snap.Result != nil {
				t.Fatal("canceled task must not expose a partial result")
			}

			// 恢复并等待完成
			if err := svc.Resume(ctx, "t"); err != nil {
				t.Fatal(err)
			}
			snap, err = svc.Wait(ctx, "t")
			if err != nil {
				t.Fatal(err)
			}
			if snap.Status != StatusCompleted {
				t.Fatalf("status after resume=%s err=%s", snap.Status, snap.Error)
			}
			if !resultsEqual(baseline, snap.Result) {
				t.Fatal("resumed result differs from uninterrupted run")
			}
		})
	}
}

// waitForCancelPoint 轮询直到任务越过指定检查点：
//   - "discover-K" ：发现第 K 个模体之后
//   - "discover"   ：整个发现阶段落盘后（progress >= 0.3）
//   - "shuffle-N"  ：第 N 次重排完成后
func waitForCancelPoint(t *testing.T, svc *Service, id, point string, shuffles, expectMotifs int) {
	t.Helper()
	target := 0.3
	if len(point) > 9 && point[:9] == "discover-" {
		n := 0
		for _, ch := range point[9:] {
			n = n*10 + int(ch-'0')
		}
		target = 0.3 * float64(n) / float64(expectMotifs)
	} else if len(point) > 8 && point[:8] == "shuffle-" {
		n := 0
		for _, ch := range point[8:] {
			n = n*10 + int(ch-'0')
		}
		target = 0.3 + 0.7*float64(n)/float64(shuffles)
	}
	for i := 0; i < 10000; i++ {
		snap, err := svc.Get(id)
		if err != nil {
			t.Fatal(err)
		}
		if snap.Progress >= target-1e-12 {
			return
		}
		if snap.Status != StatusRunning {
			t.Fatalf("task ended before cancel point: status=%s err=%s", snap.Status, snap.Error)
		}
		spinWait()
	}
	t.Fatal("timed out waiting for cancel point")
}

func TestServiceConcurrentIsolation(t *testing.T) {
	svc := newTestService(t)
	ctx := context.Background()

	motifA := []byte("ACGTATGC")
	motifB := []byte("TGCATGCA")
	seqsA, _, _ := syntheticInput(t, 11, 10, 50, motifA)
	seqsB, _, _ := syntheticInput(t, 22, 10, 50, motifB)
	cfgA := testConfig(8)
	cfgA.Shuffles = 16
	cfgB := testConfig(8)
	cfgB.Shuffles = 16

	if err := svc.Start(ctx, "A", seqsA, cfgA); err != nil {
		t.Fatal(err)
	}
	if err := svc.Start(ctx, "B", seqsB, cfgB); err != nil {
		t.Fatal(err)
	}

	snapA, err := svc.Wait(ctx, "A")
	if err != nil {
		t.Fatal(err)
	}
	snapB, err := svc.Wait(ctx, "B")
	if err != nil {
		t.Fatal(err)
	}
	if snapA.Status != StatusCompleted || snapB.Status != StatusCompleted {
		t.Fatalf("statuses A=%s B=%s (%s/%s)",
			snapA.Status, snapB.Status, snapA.Error, snapB.Error)
	}
	if snapA.Result.Motifs[0].Matrix.Consensus() != string(motifA) {
		t.Fatalf("A consensus=%q", snapA.Result.Motifs[0].Matrix.Consensus())
	}
	if snapB.Result.Motifs[0].Matrix.Consensus() != string(motifB) {
		t.Fatalf("B consensus=%q", snapB.Result.Motifs[0].Matrix.Consensus())
	}

	// 并发运行的结果必须与各自独立 Discover 完全一致
	expectA, _ := Discover(seqsA, cfgA)
	expectB, _ := Discover(seqsB, cfgB)
	if !resultsEqual(expectA, snapA.Result) {
		t.Fatal("task A result differs from standalone run")
	}
	if !resultsEqual(expectB, snapB.Result) {
		t.Fatal("task B result differs from standalone run")
	}
}

func TestServiceCrossProcessResume(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "state")
	motifBytes := []byte("GCTAGCTA")
	seqs, _, _ := syntheticInput(t, 321, 12, 50, motifBytes)
	cfg := testConfig(8)
	cfg.Shuffles = 20
	baseline, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}

	svc1, err := NewService(dir)
	if err != nil {
		t.Fatal(err)
	}
	ctx := context.Background()
	if err := svc1.Start(ctx, "cross", seqs, cfg); err != nil {
		t.Fatal(err)
	}
	waitForCancelPoint(t, svc1, "cross", "shuffle-7", cfg.Shuffles, cfg.ExpectMotifs)
	if err := svc1.Cancel("cross"); err != nil {
		t.Fatal(err)
	}
	if _, err := svc1.Wait(ctx, "cross"); err != nil {
		t.Fatal(err)
	}

	// 用一个全新的 Service 实例指向同一目录，模拟进程重启。
	svc2, err := NewService(dir)
	if err != nil {
		t.Fatal(err)
	}
	snap, err := svc2.Get("cross")
	if err != nil {
		t.Fatal(err)
	}
	if snap.Status != StatusCanceled {
		t.Fatalf("expected canceled state on disk, got %s", snap.Status)
	}
	if err := svc2.Resume(ctx, "cross"); err != nil {
		t.Fatal(err)
	}
	snap, err = svc2.Wait(ctx, "cross")
	if err != nil {
		t.Fatal(err)
	}
	if snap.Status != StatusCompleted {
		t.Fatalf("status=%s err=%s", snap.Status, snap.Error)
	}
	if !resultsEqual(baseline, snap.Result) {
		t.Fatal("cross-service resumed result differs from uninterrupted run")
	}
}
