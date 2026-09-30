package motif

import "fmt"

// Discover 执行完整的模体发现流程：
// 多轮随机初始化 EM，取最优，再做随机重排显著性检验与多重校正。
// 所有随机过程由 cfg.Seed 控制，相同输入与种子结果逐位一致。
func Discover(seqs []Seq, cfg Config) (*Result, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	if len(seqs) == 0 {
		return nil, fmt.Errorf("motif: 至少需要一条输入序列")
	}
	for i, s := range seqs {
		if len(s.Data) < cfg.MotifLength {
			return nil, fmt.Errorf("motif: 第 %d 条序列 %q 长度 %d 小于模体长度 %d",
				i, s.Name, len(s.Data), cfg.MotifLength)
		}
	}
	return DiscoverWithProgress(seqs, cfg, newCheckpoint(cfg))
}

// DiscoverWithProgress 与 Discover 相同，但通过 checkpoint 支持中断恢复。
// 恢复时不会产生部分结果或重复计算。
func DiscoverWithProgress(seqs []Seq, cfg Config, cp *Checkpoint) (*Result, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	if cp == nil {
		cp = newCheckpoint(cfg)
	}
	if err := cp.match(cfg); err != nil {
		return nil, err
	}
	if cp.Status == StatusDone && cp.result != nil {
		return cp.result, nil
	}

	bg := estimateBackground(seqs)
	filter := NewComplexityFilter(cfg.Filter)
	allowed := make([][]int, len(seqs))
	for i, s := range seqs {
		allowed[i] = filter.AllowedStarts(s.Data, cfg.MotifLength)
		if len(allowed[i]) == 0 {
			return nil, fmt.Errorf("motif: 第 %d 条序列 %q 在低复杂度过滤后没有可用位点",
				i, s.Name)
		}
	}

	motifs, err := extractMotifs(seqs, cfg, bg, allowed, cp)
	if err != nil {
		return nil, err
	}
	if err := assessSignificance(seqs, cfg, motifs, cp); err != nil {
		return nil, err
	}

	totalLen := 0
	for _, s := range seqs {
		totalLen += len(s.Data)
	}
	res := &Result{
		Seed:        cfg.Seed,
		Motifs:      motifs,
		Background:  bg,
		SequenceCnt: len(seqs),
		TotalLen:    totalLen,
	}
	cp.complete(res)
	return res, nil
}

// bestOfRounds 完成某一个模体槽位的全部 InitRounds 轮初始化（支持检查点
// 跳过已完成轮次），返回信息含量最高的模体。
func bestOfRounds(seqs []Seq, cfg Config, bg Background, allowed [][]int,
	blocked []map[int]bool, cp *Checkpoint, slot int) (*Motif, error) {

	seeds, err := seedCandidates(seqs, cfg, allowed, blocked, cp, slot)
	if err != nil {
		return nil, err
	}
	if err := cp.checkCancel(); err != nil {
		return nil, err
	}

	var best *Motif
	// 若此前已完成过部分轮次（中断恢复），先用检查点保存的槽位最优
	// 位点重建候选，保证后续新轮次能与之比较，不丢失已完成轮次的成果。
	if starts, ok := cp.bestStarts(slot); ok {
		siteBytes := make([][]byte, len(seqs))
		for i, st := range starts {
			siteBytes[i] = seqs[i].Data[st : st+cfg.MotifLength]
		}
		matrix, mErr := BuildPWM(siteBytes, bg, cfg.PWM)
		if mErr == nil {
			best = finalizeMotif(seqs, cfg, matrix, starts)
		}
	}
	for round, seed := range seeds {
		if cp.isRoundDone(slot, round) {
			continue
		}
		if err := cp.checkCancel(); err != nil {
			return nil, err
		}
		candidate := convergeFromSeed(seqs, cfg, bg, allowed, blocked, seed)
		cp.recordRound(slot, round, len(seqs), candidate)
		if err := cp.checkCancel(); err != nil {
			return nil, err
		}
		if best == nil || candidate.Information > best.Information {
			best = candidate
		}
	}
	if best == nil {
		return nil, fmt.Errorf("motif: 模体槽位 %d 未产生任何结果", slot)
	}
	return best, nil
}

// extractMotifs 从多轮最优结果中迭代抽取期望数量的模体（掩码已选位点）。
func extractMotifs(seqs []Seq, cfg Config, bg Background, allowed [][]int,
	cp *Checkpoint) ([]Motif, error) {

	blocked := make([]map[int]bool, len(seqs))
	for i := range blocked {
		blocked[i] = map[int]bool{}
	}
	if cp.Blocked != nil {
		blocked = cp.Blocked
	}
	motifs := make([]Motif, 0, cfg.ExpectedCount)
	for k := 0; k < cfg.ExpectedCount; k++ {
		if k < len(cp.Motifs) && cp.Motifs[k].Matrix != nil {
			motifs = append(motifs, cp.Motifs[k])
			applyBlocking(blocked, motifs[k], cfg.MotifLength)
			continue
		}
		best, err := bestOfRounds(seqs, cfg, bg, allowed, blocked, cp, k)
		if err != nil {
			return nil, err
		}
		best.Index = k
		motifs = append(motifs, *best)
		applyBlocking(blocked, *best, cfg.MotifLength)
		cp.recordMotif(*best, blocked)
	}
	return motifs, nil
}

func applyBlocking(blocked []map[int]bool, m Motif, motifLen int) {
	for _, occ := range m.Sites {
		for st := occ.Start - motifLen + 1; st <= occ.Start+motifLen-1; st++ {
			if st >= 0 {
				blocked[occ.SeqIndex][st] = true
			}
		}
	}
}

// estimateBackground 估计四核苷酸经验背景分布（零频率做小量平滑）。
func estimateBackground(seqs []Seq) Background {
	var counts [4]float64
	total := 0.0
	for _, s := range seqs {
		for _, c := range s.Data {
			b := baseIndex[c]
			if b >= 0 {
				counts[b]++
				total++
			}
		}
	}
	if total == 0 {
		return UniformBackground()
	}
	var bg Background
	for i := range bg {
		bg[i] = counts[i] / total
		if bg[i] == 0 {
			bg[i] = 1e-3
		}
	}
	sum := 0.0
	for _, p := range bg {
		sum += p
	}
	for i := range bg {
		bg[i] /= sum
	}
	return bg
}
