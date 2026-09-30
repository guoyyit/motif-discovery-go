package motif

import "math"

// seed 描述一个候选初始化：从某条序列上的某个窗口作为 MEME 式锚定种子。
type seed struct {
	seqIndex int
	start    int
	scout    float64
}

// seedCandidates 产生 InitRounds 个确定性候选种子。
//
// 两阶段 MEME 思路：
//  1. 侦察阶段：遍历（或分层采样）全部候选窗口，以“与全库背景的契合度”
//     的廉价代理指标打分，挑出最有希望的窗口；
//  2. 每个种子再交给 convergeFromSeed 完整收敛 EM。
//
// 随机性只用于在候选过多时的分层抽样，由 cfg.Seed 与 (槽位, 轮次)
// 确定性派生，保证可复现；相同种子无论是否中断结果一致。
func seedCandidates(seqs []Seq, cfg Config, allowed [][]int, blocked []map[int]bool,
	cp *Checkpoint, slot int) ([]seed, error) {
	rng := StreamFromSeed(cfg.Seed, seedStreamID(slot))

	// k-mer 出现集合：kmerPresence[word] 为包含该 word 的序列编号集合。
	kmerPresence := map[uint64]map[int]struct{}{}
	for i, s := range seqs {
		for st := 0; st+cfg.MotifLength <= len(s.Data); st++ {
			word := encodeKmer(s.Data[st : st+cfg.MotifLength])
			set, ok := kmerPresence[word]
			if !ok {
				set = map[int]struct{}{}
				kmerPresence[word] = set
			}
			set[i] = struct{}{}
		}
	}

	type candidate = scoredSeed
	var candidates []candidate

	// 侦察指标：该精确窗口在其余序列中的出现序列数（富集潜力的廉价代理）。
	for i, s := range seqs {
		starts := availableStarts(allowed[i], blocked[i])
		if cp != nil {
			if err := cp.checkCancel(); err != nil {
				return nil, err
			}
		}
		for _, st := range starts {
			w := s.Data[st : st+cfg.MotifLength]
			word := encodeKmer(w)
			score := float64(len(kmerPresence[word]) - 1)
			candidates = append(candidates, candidate{
				seed:  seed{seqIndex: i, start: st},
				score: score,
			})
		}
	}

	// 候选过多时分两层处理：保留每个侦察分桶的代表，并辅以确定性随机抽样，
	// 保证 InitRounds 个种子既有高潜力候选，也覆盖搜索空间。
	k := cfg.InitRounds
	if len(candidates) <= k {
		out := make([]seed, len(candidates))
		for i, c := range candidates {
			out[i] = c.seed
		}
		return out, nil
	}

	topK := topSeeds(candidates, k/2+1)
	out := make([]seed, 0, k)
	out = append(out, topK...)
	for len(out) < k {
		idx := rng.Intn(len(candidates))
		out = append(out, candidates[idx].seed)
	}
	return out, nil
}

// encodeKmer 将长度不超过 32 的 ACGT k-mer 编码为 2bit 整数。
func encodeKmer(w []byte) uint64 {
	var x uint64
	for _, c := range w {
		x <<= 2
		switch c {
		case 'C':
			x |= 1
		case 'G':
			x |= 2
		case 'T':
			x |= 3
		}
	}
	return x
}

type scoredSeed struct {
	seed  seed
	score float64
}

func topSeeds(items []scoredSeed, n int) []seed {
	for i := 0; i < n && i < len(items); i++ {
		maxIdx := i
		for j := i + 1; j < len(items); j++ {
			if items[j].score > items[maxIdx].score {
				maxIdx = j
			}
		}
		items[i], items[maxIdx] = items[maxIdx], items[i]
	}
	out := make([]seed, 0, n)
	for i := 0; i < n && i < len(items); i++ {
		if items[i].score <= 0 {
			break
		}
		out = append(out, items[i].seed)
	}
	return out
}

// convergeFromSeed 从单个 MEME 锚定种子出发执行硬分配 EM 直到收敛。
// 初始矩阵由锚定窗口主导（8 份）加每条序列一个随机位点（各 1 份），
// 随后交替执行“每序列取最佳位点”（E 步）与“由位点重建 PWM”（M 步）。
func convergeFromSeed(seqs []Seq, cfg Config, bg Background, allowed [][]int,
	blocked []map[int]bool, sd seed) *Motif {

	rng := StreamFromSeed(cfg.Seed, emJitterStreamID(sd.seqIndex, sd.start))

	current := make([]int, len(seqs))
	initSites := make([][]byte, 0, len(seqs)+8)
	for rep := 0; rep < 8; rep++ {
		initSites = append(initSites,
			seqs[sd.seqIndex].Data[sd.start:sd.start+cfg.MotifLength])
	}
	for i := range seqs {
		starts := availableStarts(allowed[i], blocked[i])
		pos := starts[rng.Intn(len(starts))]
		if i == sd.seqIndex {
			pos = sd.start
		}
		current[i] = pos
		initSites = append(initSites, seqs[i].Data[pos:pos+cfg.MotifLength])
	}

	matrix, err := BuildPWM(initSites, bg, cfg.PWM)
	if err != nil {
		// 种子始终是合法窗口，这里理论上不会触发。
		panic(err)
	}

	prevMean := math.Inf(-1)
	for iter := 0; iter < cfg.MaxEMIterations; iter++ {
		next := make([]int, len(seqs))
		meanNorm := 0.0
		for i, s := range seqs {
			starts := availableStarts(allowed[i], blocked[i])
			bestStart := current[i]
			bestNorm := math.Inf(-1)
			for _, st := range starts {
				ns := matrix.NormalizedScore(s.Data[st : st+cfg.MotifLength])
				if ns > bestNorm {
					bestNorm = ns
					bestStart = st
				}
			}
			next[i] = bestStart
			meanNorm += bestNorm
		}
		meanNorm /= float64(len(seqs))

		siteBytes := make([][]byte, len(seqs))
		for i, st := range next {
			siteBytes[i] = seqs[i].Data[st : st+cfg.MotifLength]
		}
		matrix, err = BuildPWM(siteBytes, bg, cfg.PWM)
		if err != nil {
			panic(err)
		}
		stable := sameAssignments(current, next)
		current = next
		if stable {
			break
		}
		if math.Abs(meanNorm-prevMean) <= cfg.ConvergenceTol*math.Max(1, math.Abs(prevMean)) {
			break
		}
		prevMean = meanNorm
	}

	return finalizeMotif(seqs, cfg, matrix, current)
}

// finalizeMotif 由最终 PWM 与每条序列的最佳位点构造 Motif。
func finalizeMotif(seqs []Seq, cfg Config, matrix *PWM, starts []int) *Motif {
	occs := make([]SiteOccurrence, len(seqs))
	mean := 0.0
	enriched := 0
	info := 0.0
	for i, st := range starts {
		ns := matrix.NormalizedScore(seqs[i].Data[st : st+cfg.MotifLength])
		occs[i] = SiteOccurrence{SeqIndex: i, Start: st, Score: ns}
		mean += ns
		if ns >= cfg.PWM.Threshold {
			enriched++
		}
	}
	mean /= float64(len(seqs))
	for _, v := range matrix.Info {
		info += v
	}
	return &Motif{
		Matrix:      matrix,
		Consensus:   matrix.Consensus,
		Score:       mean,
		Enrichment:  float64(enriched) / float64(len(seqs)),
		Sites:       occs,
		Information: info,
	}
}

func availableStarts(allowed []int, blocked map[int]bool) []int {
	if len(blocked) == 0 {
		return allowed
	}
	out := make([]int, 0, len(allowed))
	for _, st := range allowed {
		if !blocked[st] {
			out = append(out, st)
		}
	}
	return out
}

func sameAssignments(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func seedStreamID(slot int) uint64 {
	return 4_000_000_007 + uint64(slot)*1_000_003
}

func emJitterStreamID(seqIndex, start int) uint64 {
	return 5_000_000_009 + uint64(seqIndex)*7919 + uint64(start)*31
}
