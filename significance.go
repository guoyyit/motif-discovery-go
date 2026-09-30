package motif

import (
	"sort"
)

// SigConfig 控制基于随机重排零分布的显著性检验。
type SigConfig struct {
	// Permutations 为零分布随机重排次数 B。
	Permutations int
	// Alpha 为显著性阈值（校正后）。
	Alpha float64
	// Method 为多重检验校正方法，目前支持 "BH"（Benjamini-Hochberg）。
	Method string
}

// DefaultSigConfig 返回默认的显著性检验配置。
func DefaultSigConfig() SigConfig {
	return SigConfig{
		Permutations: 100,
		Alpha:        0.05,
		Method:       "BH",
	}
}

// NullDistribution 保存随机重排零分布的统计量样本。
type NullDistribution struct {
	Stats []float64
}

// motifNullStatistic 以发现到的 PWM 为检验统计量：对每条（已随机重排的）
// 序列，在通过低复杂度过滤的窗口中取最佳归一化得分，再对所有序列取均值。
// 观测值高而重排值低，说明该模体的富集无法由序列碱基组成解释。
func motifNullStatistic(matrix *PWM, shuffled [][]byte, allowed [][]int) float64 {
	total := 0.0
	used := 0
	for i, data := range shuffled {
		best := -1.0
		for _, st := range allowed[i] {
			if st+matrix.Length > len(data) {
				continue
			}
			ns := matrix.NormalizedScore(data[st : st+matrix.Length])
			if ns > best {
				best = ns
			}
		}
		if best >= 0 {
			total += best
			used++
		}
	}
	if used == 0 {
		return 0
	}
	return total / float64(used)
}

// assessSignificance 为每个已发现模体构建随机重排零分布。
//
// 统计语义（关键）：每一次重排都重新执行完整的多轮初始化发现流程，
// 因此零分布自动包含了“从数据中挑选最优模体”这一选择过程（含赢家诅咒），
// 经验 p 值在 H0 下近似均匀，而不会因事后选择产生反保守偏差。
// 重排保持每条序列的碱基组成（Fisher-Pitman 式置换）。
//
// 对期望 K 个模体，使用顺序检验：第 k 次重排发现的第 k 个模体得分
// 作为观测第 k 个模体的零分布，最后统一做 BH FDR 校正。
//
// 整个过程支持检查点恢复：每个 (槽位, 重排编号) 的随机流由种子与编号
// 确定性派生，已完成的重放不会重复计算；结果仅在全部完成后一次性写回。
func assessSignificance(seqs []Seq, cfg Config, motifs []Motif, cp *Checkpoint) error {
	filter := NewComplexityFilter(cfg.Filter)
	for k := range motifs {
		stats := cp.nullStats(k)
		for b := len(stats); b < cfg.Significance.Permutations; b++ {
			if err := cp.checkCancel(); err != nil {
				return err
			}
			rng := StreamFromSeed(cfg.Seed, nullStreamID(k, b))
			permSeqs := make([]Seq, len(seqs))
			allowed := make([][]int, len(seqs))
			for i, s := range seqs {
				permSeqs[i] = Seq{Name: s.Name, Desc: s.Desc, Data: rng.Shuffled(s.Data)}
				allowed[i] = filter.AllowedStarts(permSeqs[i].Data, cfg.MotifLength)
				if len(allowed[i]) == 0 {
					allowed[i] = []int{0}
					if len(permSeqs[i].Data) > cfg.MotifLength {
						allowed[i][0] = rng.Intn(len(permSeqs[i].Data) - cfg.MotifLength + 1)
					}
				}
			}
			bg := estimateBackground(permSeqs)
			nullMotifs, err := rediscoverWithoutTest(permSeqs, cfg, bg, allowed, k)
			if err != nil {
				return err
			}
			stats = append(stats, nullMotifs[k].Information)
			cp.recordNull(k, stats)
		}

		motifs[k].PValue = EmpiricalPValue(motifs[k].Information, stats)
		cp.recordMotifUpdate(k, motifs[k])
	}

	qvalues := BHAdjust(pvaluesOf(motifs))
	for k := range motifs {
		motifs[k].AdjustedP = qvalues[k]
		motifs[k].Significant = qvalues[k] <= cfg.Significance.Alpha
		cp.recordMotifUpdate(k, motifs[k])
	}
	return nil
}

// rediscoverWithoutTest 在重排数据上运行发现的前 (slot+1) 个槽位，
// 跳过显著性检验与检查点记账（零分布内部调用，必须无副作用、确定）。
func rediscoverWithoutTest(seqs []Seq, cfg Config, bg Background, allowed [][]int, slot int) ([]Motif, error) {
	blocked := make([]map[int]bool, len(seqs))
	for i := range blocked {
		blocked[i] = map[int]bool{}
	}
	out := make([]Motif, 0, slot+1)
	for k := 0; k <= slot; k++ {
		var best *Motif
		seeds, err := seedCandidates(seqs, cfg, allowed, blocked, nil, k)
		if err != nil {
			return nil, err
		}
		for _, sd := range seeds {
			candidate := convergeFromSeed(seqs, cfg, bg, allowed, blocked, sd)
			if best == nil || candidate.Information > best.Information {
				best = candidate
			}
		}
		best.Index = k
		out = append(out, *best)
		applyBlocking(blocked, *best, cfg.MotifLength)
	}
	return out, nil
}

// PermutedNull 对序列集合做组成保留随机重排并在重排数据上重新运行完整
// 发现流程（一轮），返回指定槽位的最优模体得分。公开 API，便于外部
// 自定义/验证零分布。
func PermutedNull(seqs []Seq, cfg Config, rng *RNG, slot int) (float64, error) {
	filter := NewComplexityFilter(cfg.Filter)
	permSeqs := make([]Seq, len(seqs))
	allowed := make([][]int, len(seqs))
	for i, s := range seqs {
		permSeqs[i] = Seq{Name: s.Name, Data: rng.Shuffled(s.Data)}
		allowed[i] = filter.AllowedStarts(permSeqs[i].Data, cfg.MotifLength)
	}
	bg := estimateBackground(permSeqs)
	ms, err := rediscoverWithoutTest(permSeqs, cfg, bg, allowed, slot)
	if err != nil {
		return 0, err
	}
	return ms[slot].Information, nil
}

// EmpiricalPValue 基于零分布样本计算经验 p 值
// p = (1 + #{s >= obs}) / (1 + B)，含 +1 校正，故永远不会得到 0。
func EmpiricalPValue(observed float64, null []float64) float64 {
	ge := 0
	for _, s := range null {
		if s >= observed {
			ge++
		}
	}
	return float64(1+ge) / float64(1+len(null))
}

// BHAdjust 对 p 值做 Benjamini-Hochberg FDR 多重检验校正，
// 返回与输入等长、等序的校正后 q 值。
func BHAdjust(pvalues []float64) []float64 {
	n := len(pvalues)
	q := make([]float64, n)
	if n == 0 {
		return q
	}
	order := make([]int, n)
	for i := range order {
		order[i] = i
	}
	sort.Slice(order, func(a, b int) bool {
		return pvalues[order[a]] < pvalues[order[b]]
	})
	running := 1.0
	for rank := n; rank >= 1; rank-- {
		idx := order[rank-1]
		v := pvalues[idx] * float64(n) / float64(rank)
		if v > running {
			v = running
		}
		if v > 1 {
			v = 1
		}
		running = v
		q[idx] = v
	}
	return q
}

func pvaluesOf(motifs []Motif) []float64 {
	out := make([]float64, len(motifs))
	for i := range motifs {
		out[i] = motifs[i].PValue
	}
	return out
}

// nullStreamID 为 (模体槽位, 重排编号) 派生确定性随机流。
// 大基数偏移避免与侦察/EM 阶段的随机流空间重叠。
func nullStreamID(slot, perm int) uint64 {
	return 1_000_000_007 + uint64(slot)*1_000_003 + uint64(perm)
}
