package motif

import "testing"

// TestSignificanceCalibration 在无模体的随机序列集合上多次重复实验，
// 经验假阳性率不应显著超过名义水平 alpha。
// 使用较小的重排次数（B=39）使最小可达 p 值为 1/40=0.025，
// 在 alpha=0.05 下若零假设成立，显著比例应处于较低水平。
func TestSignificanceCalibration(t *testing.T) {
	const trials = 12
	const alpha = 0.05
	cfg := testConfig(8)
	cfg.Shuffles = 39
	cfg.SignificanceAlpha = alpha

	falsePositives := 0
	anySignificant := 0
	for trial := 0; trial < trials; trial++ {
		seqs, _, _ := syntheticInput(t, uint64(700+trial), 12, 50, nil)
		result, err := Discover(seqs, cfg)
		if err != nil {
			t.Fatal(err)
		}
		sigInTrial := false
		for _, m := range result.Motifs {
			if m.PValue <= 0 || m.PValue > 1 {
				t.Fatalf("p-value out of range: %g", m.PValue)
			}
			if m.AdjustedPValue < m.PValue-1e-12 {
				t.Fatalf("adjusted p (%g) smaller than raw p (%g)",
					m.AdjustedPValue, m.PValue)
			}
			if m.Significant {
				sigInTrial = true
				falsePositives++
			}
		}
		if sigInTrial {
			anySignificant++
		}
	}

	// 以试验为单位的假阳性率：12 次试验中至多允许 3 次出现显著模体。
	// 该界比名义 5% 宽松，是为了在极小样本下避免测试抖动，
	// 同时仍能捕捉严重的校准失效（如任何试验都显著）。
	if anySignificant > 3 {
		t.Fatalf("calibration failure: %d/%d null trials produced significant motifs",
			anySignificant, trials)
	}
	t.Logf("null trials with a significant motif: %d/%d (motif-level FPs=%d)",
		anySignificant, trials, falsePositives)
}

// TestSignificanceStrongSignal 验证强植入信号几乎总是显著。
func TestSignificanceStrongSignal(t *testing.T) {
	cfg := testConfig(8)
	cfg.Shuffles = 39
	passed := 0
	const trials = 5
	for trial := 0; trial < trials; trial++ {
		seqs, _, _ := syntheticInput(t, uint64(800+trial), 15, 50, []byte("ACGTATGC"))
		result, err := Discover(seqs, cfg)
		if err != nil {
			t.Fatal(err)
		}
		if result.Motifs[0].Significant {
			passed++
		}
	}
	if passed < 4 {
		t.Fatalf("strong signal detected in only %d/%d trials", passed, trials)
	}
}
