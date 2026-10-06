package motif

import "testing"

// runNullWithRounds 用固定的观测拟合结果，按指定初始化轮数运行零分布拟合，
// 返回累积的零分布统计量。用于隔离观测/零分布两侧的拟合强度。
func runNullWithRounds(t *testing.T, seqs []Sequence, cfg Config, rounds int) [][]float64 {
	t.Helper()
	state, err := newDiscoveryState(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := state.runDiscover(nil); err != nil {
		t.Fatal(err)
	}
	// 观测拟合完成后单独调整零分布拟合强度。
	state.cfg.InitRounds = rounds
	if err := state.runShuffles(nil); err != nil {
		t.Fatal(err)
	}
	return state.nullScores
}

// TestNullFitUsesConfiguredStrength 零分布拟合必须使用与观测相同的
// 拟合强度（cfg.InitRounds），而不是独立的更弱强度：
// 改变 InitRounds 必须改变零分布；相同 InitRounds 必须逐位一致。
func TestNullFitUsesConfiguredStrength(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 909, 10, 50, []byte("ACGTATGC"))
	cfg := testConfig(8)
	cfg.Shuffles = 8

	null1 := runNullWithRounds(t, seqs, cfg, 1)
	null6 := runNullWithRounds(t, seqs, cfg, 6)

	// 若零分布忽略 cfg.InitRounds（例如硬编码为固定轮数），
	// 两种设置下的零分布会完全相同。
	same := len(null1) == len(null6)
	if same {
	outer:
		for k := range null1 {
			if len(null1[k]) != len(null6[k]) {
				same = false
				break
			}
			for i := range null1[k] {
				if null1[k][i] != null6[k][i] {
					same = false
					break outer
				}
			}
		}
	}
	if same {
		t.Fatal("null distribution ignores cfg.InitRounds: fitting strength is not shared with observed fit")
	}

	// 相同强度重复运行必须逐位一致。
	null1again := runNullWithRounds(t, seqs, cfg, 1)
	if len(null1) != len(null1again) {
		t.Fatal("null distribution size changed between identical runs")
	}
	for k := range null1 {
		for i := range null1[k] {
			if null1[k][i] != null1again[k][i] {
				t.Fatalf("null score mismatch at motif %d shuffle %d: %g vs %g",
					k, i, null1[k][i], null1again[k][i])
			}
		}
	}
}

// TestDiscoverBitwiseReproducible 相同输入与相同种子下重复运行必须
// 逐位一致，包括 PWM 矩阵的全部浮点系数（比 resultsEqual 更严格）。
func TestDiscoverBitwiseReproducible(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 910, 12, 50, []byte("GCTAGCTA"))
	cfg := testConfig(8)
	cfg.ExpectMotifs = 2
	cfg.Shuffles = 20

	r1, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	r2, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if !resultsEqual(r1, r2) {
		t.Fatal("repeated runs differ")
	}
	for k := range r1.Motifs {
		a, b := r1.Motifs[k].Matrix, r2.Motifs[k].Matrix
		for pos := range a.LogOdds {
			for base := range a.LogOdds[pos] {
				if a.LogOdds[pos][base] != b.LogOdds[pos][base] {
					t.Fatalf("motif %d log-odds differ at (%d,%d): %g vs %g",
						k, pos, base, a.LogOdds[pos][base], b.LogOdds[pos][base])
				}
			}
		}
	}
}
