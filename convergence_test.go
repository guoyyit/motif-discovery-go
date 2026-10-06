package motif

import (
	"errors"
	"strings"
	"testing"
	"time"
)

// within 在超时前完成 f 才算通过，用于“有限步内结束”的守护。
func within(t *testing.T, d time.Duration, what string, f func()) {
	t.Helper()
	done := make(chan struct{})
	go func() {
		defer close(done)
		f()
	}()
	select {
	case <-done:
	case <-time.After(d):
		t.Fatalf("%s did not finish within %v", what, d)
	}
}

// TestEMDegenerateInputTerminates 退化输入（全部窗口完全相同）下
// EM 必须在有限步内收敛，不得挂起。
func TestEMDegenerateInputTerminates(t *testing.T) {
	seqs := [][]byte{
		[]byte(strings.Repeat("ACGT", 8)),
		[]byte(strings.Repeat("ACGT", 8)),
		[]byte(strings.Repeat("ACGT", 8)),
	}
	windows := enumerateWindows(seqs, 8, nil)
	within(t, 10*time.Second, "EM on degenerate input", func() {
		model, err := fitOneMotif(windows, 8, DefaultBackground(), 0.5, 1e-4, 5, 500,
			7, globalStream(300, 0))
		if err != nil {
			t.Fatalf("degenerate input should converge, got %v", err)
		}
		if model == nil || model.matrix == nil {
			t.Fatal("converged fit returned no model")
		}
	})
}

// TestEMMaxIterationsNotConverged 迭代上限为 1 时任何拟合都不可能收敛，
// 必须以可区分的 ErrNotConverged 返回，且不得返回半成品模型。
func TestEMMaxIterationsNotConverged(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 77, 8, 40, nil)
	raw := make([][]byte, len(seqs))
	for i, s := range seqs {
		raw[i] = s.Data
	}
	windows := enumerateWindows(raw, 8, nil)
	within(t, 10*time.Second, "EM with MaxIterations=1", func() {
		model, err := fitOneMotif(windows, 8, DefaultBackground(), 0.5, 1e-4, 5, 1,
			7, globalStream(301, 0))
		if !errors.Is(err, ErrNotConverged) {
			t.Fatalf("expected ErrNotConverged, got model=%v err=%v", model, err)
		}
		if model != nil {
			t.Fatal("non-converged fit must not return a partial model")
		}
	})
}

// TestEMExtremeToleranceTerminates 极端收敛容差（远超浮点精度）下
// 拟合必须在迭代上限内结束：要么收敛，要么明确不收敛，不得挂起。
func TestEMExtremeToleranceTerminates(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 78, 8, 40, []byte("ACGTATGC"))
	raw := make([][]byte, len(seqs))
	for i, s := range seqs {
		raw[i] = s.Data
	}
	windows := enumerateWindows(raw, 8, nil)
	within(t, 30*time.Second, "EM with tol=1e-300", func() {
		_, err := fitOneMotif(windows, 8, DefaultBackground(), 0.5, 1e-300, 3, 200,
			7, globalStream(302, 0))
		if err != nil && !errors.Is(err, ErrNotConverged) {
			t.Fatalf("unexpected error: %v", err)
		}
	})
}

// TestDiscoverReportsNotConverged Discover 层面：迭代上限过小时
// 任务必须以可区分的不收敛错误结束，而不是挂起或返回部分结果。
func TestDiscoverReportsNotConverged(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 79, 8, 40, []byte("ACGTATGC"))
	cfg := testConfig(8)
	cfg.Shuffles = 2
	cfg.MaxIterations = 1
	within(t, 30*time.Second, "Discover with MaxIterations=1", func() {
		result, err := Discover(seqs, cfg)
		if !errors.Is(err, ErrNotConverged) {
			t.Fatalf("expected ErrNotConverged, got result=%v err=%v", result, err)
		}
		if result != nil {
			t.Fatal("failed discovery must not return a partial result")
		}
	})
}

// TestEMConvergenceConfigurable 收敛容差与迭代上限均可配置：
// 宽松的容差应当让拟合正常收敛并给出完整模型。
func TestEMConvergenceConfigurable(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 80, 10, 40, []byte("ACGTATGC"))
	raw := make([][]byte, len(seqs))
	for i, s := range seqs {
		raw[i] = s.Data
	}
	windows := enumerateWindows(raw, 8, nil)
	model, err := fitOneMotif(windows, 8, DefaultBackground(), 0.5, 0.5, 5, 50,
		7, globalStream(303, 0))
	if err != nil {
		t.Fatalf("loose tolerance should converge: %v", err)
	}
	if model.matrix.Width != 8 {
		t.Fatalf("unexpected model: %+v", model)
	}
}
