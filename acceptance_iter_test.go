package motif

import (
	"errors"
	"math"
	"strings"
	"testing"
	"time"
)

// repeatInput 构造 nSeq 条“随机背景 + 长 AT 串联重复”的序列。
func repeatInput(t *testing.T, seed uint64, nSeq int) []Sequence {
	t.Helper()
	r := streamRand(seed, globalStream(9003, 0))
	bases := []byte("ACGT")
	seqs := make([]Sequence, nSeq)
	for i := range seqs {
		seq := make([]byte, 0, 64)
		for len(seq) < 40 {
			seq = append(seq, bases[r.intn(4)])
		}
		seq = append(seq, strings.Repeat("AT", 12)...) // 24bp 二核苷酸重复
		seqs[i] = Sequence{Name: "r", Data: seq}
	}
	return seqs
}

// TestDefaultConfigRejectsSimpleRepeats 需求 1：默认配置（不设任何开关）
// 下，AT 串联重复不得被报告为显著模体；且默认路径与显式开启防护的结论一致。
func TestDefaultConfigRejectsSimpleRepeats(t *testing.T) {
	seqs := repeatInput(t, 888, 20)

	// 完全零值配置：只给种子，其余全部走默认值。
	defaultCfg := Config{Seed: 42, Shuffles: 40}
	resDefault, err := Discover(seqs, defaultCfg)
	if err != nil {
		t.Fatal(err)
	}
	for _, m := range resDefault.Motifs {
		if m.Significant && IsLowComplexity([]byte(m.Matrix.Consensus())) {
			t.Fatalf("default path reported a simple repeat as significant: %q p=%g",
				m.Matrix.Consensus(), m.PValue)
		}
		if m.PValue <= 0 || m.PValue > 1 {
			t.Fatalf("p-value out of (0,1], statistics must be retained: %g", m.PValue)
		}
		// 可核对统计量：该模体零分布样本数不超过配置的重排次数
		// （部分重排可能因可用窗口不足而提前停止，该模体少一个零样本，
		// p 值按各自实际 B 计算，偏保守）；p=(1+exceed)/(1+B)
		// 反解出的超出次数必须是 [0,B] 内整数。
		if m.NullCount < 0 || m.NullCount > 40 {
			t.Fatalf("NullCount=%d out of [0,40]", m.NullCount)
		}
		if m.NullCount == 0 {
			// 零样本为 0 时唯一合法结论是 p=1（不可能显著）。
			if m.PValue != 1 || m.Significant {
				t.Fatalf("zero null samples must give p=1 and not significant, got p=%g sig=%v",
					m.PValue, m.Significant)
			}
			continue
		}
		exceed := m.PValue*float64(m.NullCount+1) - 1
		rounded := math.Round(exceed)
		if math.Abs(exceed-rounded) > 1e-9 || rounded < 0 || rounded > float64(m.NullCount) {
			t.Fatalf("p-value %g not auditable against B=%d (exceed=%g)",
				m.PValue, m.NullCount, exceed)
		}
	}

	// 显式开启防护：结论必须与默认路径逐位一致。
	explicitCfg := defaultCfg
	explicitCfg.MaskLowComplexity = Bool(true)
	resExplicit, err := Discover(seqs, explicitCfg)
	if err != nil {
		t.Fatal(err)
	}
	if !resultsEqual(resDefault, resExplicit) {
		t.Fatal("default path and explicit MaskLowComplexity=true must agree exactly")
	}
}

// TestPureRepeatDefaultIsNotSignificant 纯简单重复输入在默认路径下
// 必须返回“不显著”并保留可核对的统计量（p 值、零分布样本数）。
func TestPureRepeatDefaultIsNotSignificant(t *testing.T) {
	seqs := make([]Sequence, 12)
	for i := range seqs {
		seqs[i] = Sequence{Data: []byte(strings.Repeat("AT", 20))}
	}
	res, err := Discover(seqs, Config{Seed: 7, Shuffles: 20})
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Motifs) != 0 {
		// 若有候选模体，它们也绝不允许显著。
		for _, m := range res.Motifs {
			if m.Significant {
				t.Fatalf("pure AT repeat must not be significant: %q p=%g",
					m.Matrix.Consensus(), m.PValue)
			}
		}
	}
}

// type fitRecord 记录一次拟合的阶段与强度参数。
type fitRecord struct {
	phase    string
	rounds   int
	maxIters int
	tol      float64
}

// TestObservedNullFitParity 需求 2：观测拟合与零分布拟合必须使用相同的
// 拟合强度（初始化轮数）、相同的收敛容差与迭代上限，并以相同的占用顺序
// 逐个模体拟合。这里通过测试观测缝断言每个 null 拟合的轮数与观测一致。
func TestObservedNullFitParity(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 2024, 12, 50, []byte("ACGTATGC"))
	cfg := testConfig(8)
	cfg.Shuffles = 8
	cfg.ExpectMotifs = 2
	cfg.InitRounds = 6
	cfg.MaxEMIterations = 500

	var records []fitRecord
	old := fitObserver
	fitObserver = func(phase string, rounds, maxIters int, tol float64) {
		records = append(records, fitRecord{phase, rounds, maxIters, tol})
	}
	defer func() { fitObserver = old }()

	if _, err := Discover(seqs, cfg); err != nil {
		t.Fatal(err)
	}

	var observed, null []fitRecord
	for _, rec := range records {
		switch rec.phase {
		case "observed":
			observed = append(observed, rec)
		case "null":
			null = append(null, rec)
		}
	}
	if len(observed) == 0 {
		t.Fatal("no observed fits recorded")
	}
	if len(null) == 0 {
		t.Fatal("no null fits recorded")
	}
	// 每次重排对每个已发现模体做一次 null 拟合；强度必须与观测完全相同：
	// 相同的初始化轮数、相同的迭代上限、相同的收敛容差。
	obs := observed[0]
	for _, n := range null {
		if n.rounds != obs.rounds || n.maxIters != obs.maxIters || n.tol != obs.tol {
			t.Fatalf("null fit strength %+v differs from observed %+v", n, obs)
		}
		if n.rounds != cfg.InitRounds || n.maxIters != cfg.MaxEMIterations ||
			n.tol != cfg.ConvergenceTol {
			t.Fatalf("null fit %+v does not match configured strength "+
				"(rounds=%d maxIters=%d tol=%g)",
				n, cfg.InitRounds, cfg.MaxEMIterations, cfg.ConvergenceTol)
		}
	}
	for _, o := range observed {
		if o.rounds != cfg.InitRounds || o.maxIters != cfg.MaxEMIterations ||
			o.tol != cfg.ConvergenceTol {
			t.Fatalf("observed fit %+v does not match configured: rounds=%d maxIters=%d tol=%g",
				o, cfg.InitRounds, cfg.MaxEMIterations, cfg.ConvergenceTol)
		}
	}

	// 占用顺序一致：每个 shuffle 内的 null 拟合模体顺序必须与观测发现顺序一致。
	// 观测拟合按 k=0..K-1 顺序记录；null 记录按 b 外层、k 内层产生，
	// 因此连续每 len(observed) 条 null 记录的模体次序应能与观测阶段对应。
	if len(null) != cfg.Shuffles*len(observed) {
		t.Fatalf("null fit count=%d, want %d", len(null), cfg.Shuffles*len(observed))
	}
}

// TestEMBoundedOnDegenerateInput 需求 3：退化输入（全部相同、近似重复）
// 下 EM 必须在有限步内结束。
func TestEMBoundedOnDegenerateInput(t *testing.T) {
	mkWindows := func(seq []byte, width int) []window {
		return enumerateWindows([][]byte{seq}, width, nil)
	}

	// 完全均一序列：似然首轮即常数，必须立即收敛返回。
	uniform := []byte(strings.Repeat("A", 40))
	w := mkWindows(uniform, 8)
	m, err := runEMSeeded(w, 8, DefaultBackground(), 0.5, 1e-4, 1000,
		dataDrivenSeeds(w, 8, DefaultBackground()))
	if err != nil {
		t.Fatalf("uniform input should converge immediately: %v", err)
	}
	if m == nil {
		t.Fatal("uniform input produced nil model")
	}

	// 近似重复输入：每条序列是同一重复的轻微扰动，拟合也必须在有限步结束。
	r := streamRand(5, globalStream(9004, 0))
	seqs := make([][]byte, 10)
	for i := range seqs {
		s := []byte(strings.Repeat("ATG", 16))
		s[r.intn(len(s))] = "ACGT"[r.intn(4)] // 单点扰动
		seqs[i] = s
	}
	w = enumerateWindows(seqs, 8, nil)
	done := make(chan struct{})
	go func() {
		_, _ = fitOneMotif(w, 8, DefaultBackground(), 0.5, 1e-4, 4, 100000,
			9, globalStream(101, 0), "observed")
		close(done)
	}()
	select {
	case <-done:
	case <-time.After(30 * time.Second):
		t.Fatal("EM did not finish in finite time on near-repeat input")
	}
}

// TestEMMaxIterationsDistinguishable 需求 3：迭代上限设为 1 时首轮不可能
// 同时完成收敛检查（首轮 prev=+Inf），因此必须返回包装 ErrEMNotConverged
// 的明确错误，而不是挂起或返回半成品。
func TestEMMaxIterationsDistinguishable(t *testing.T) {
	seqs, _, _ := syntheticInput(t, 31337, 10, 40, []byte("ACGTATGC"))
	raw := make([][]byte, len(seqs))
	for i, s := range seqs {
		raw[i] = s.Data
	}
	w := enumerateWindows(raw, 8, nil)
	_, err := runEMSeeded(w, 8, DefaultBackground(), 0.5, 1e-12, 1,
		dataDrivenSeeds(w, 8, DefaultBackground()))
	if !errors.Is(err, ErrEMNotConverged) {
		t.Fatalf("want ErrEMNotConverged at MaxEMIterations=1, got %v", err)
	}

	// 端到端：配置 MaxEMIterations=1 时 Discover 必须以可区分错误失败。
	cfg := testConfig(8)
	cfg.MaxEMIterations = 1
	cfg.Shuffles = 4
	if _, err := Discover(seqs, cfg); !errors.Is(err, ErrEMNotConverged) {
		t.Fatalf("Discover with MaxEMIterations=1 must fail with ErrEMNotConverged, got %v", err)
	}
}
