package motif

import (
	"errors"
	"fmt"
	"math"
	"sort"
)

// ErrEMNotConverged 表示 EM 在配置的迭代上限内未达到收敛条件。
// 调用方可用 errors.Is 区分“明确的不收敛结论”与其他错误：
// 达到上限的拟合轮次会被丢弃；若一次模体拟合的全部初始化轮次
// （含确定性种子轮）都未收敛，fitOneMotif 返回包装该哨兵的错误，
// 不会返回半成品模型，也不会无限挂起。
var ErrEMNotConverged = errors.New("motif: EM did not converge within MaxEMIterations")

// window 表示一条序列中一个可打分的有效片段。
type window struct {
	seqIndex int
	offset   int
	data     []byte
}

// emModel 是一次 EM 拟合的中间结果。
type emModel struct {
	matrix *PWM
	sites  []Site
	loglik float64
}

// fitObserver 是仅供测试使用的观测缝：每次 fitOneMotif 调用时记录
// 拟合阶段（"observed"/"null"）、实际使用的初始化轮数、迭代上限与
// 收敛容差，用于断言观测拟合与零分布拟合的强度完全一致。生产代码中保持为 nil。
var fitObserver func(phase string, rounds, maxIters int, tol float64)

func reportFit(phase string, rounds, maxIters int, tol float64) {
	if fitObserver != nil {
		fitObserver(phase, rounds, maxIters, tol)
	}
}

// enumerateWindows 枚举所有长度为 width、不含未知碱基 N 的片段。
// blocked 为可选的按序列索引的已占用区间集合，落入其中的片段会被跳过。
func enumerateWindows(seqs [][]byte, width int, blocked map[int][]interval) []window {
	var windows []window
	for si, seq := range seqs {
		for offset := 0; offset+width <= len(seq); offset++ {
			valid := true
			for j := 0; j < width; j++ {
				if seq[offset+j] == 'N' {
					valid = false
					break
				}
			}
			if !valid {
				continue
			}
			if overlapsAny(blocked[si], offset, offset+width) {
				continue
			}
			windows = append(windows, window{
				seqIndex: si,
				offset:   offset,
				data:     seq[offset : offset+width],
			})
		}
	}
	return windows
}

type interval struct {
	start int
	end   int
}

func overlapsAny(ivs []interval, start, end int) bool {
	for _, iv := range ivs {
		if start < iv.end && iv.start < end {
			return true
		}
	}
	return false
}

// fitOneMotif 对给定窗口集合运行 rounds 次随机初始化的 EM，返回似然最高的模型。
// 每轮使用独立的确定性子流（baseStream + round），逐轮顺序与中断恢复无关。
//
// 观测拟合（phase="observed"）与零分布拟合（phase="null"）必须以完全相同的
// 强度调用本函数：相同的 rounds、tol 与 maxIters。零分布不得使用更弱的拟合，
// 否则显著性结论会随拟合强度漂移。
func fitOneMotif(windows []window, width int, bg Background, pseudocount, tol float64,
	rounds, maxIters int, seed uint64, baseStream uint64, phase string) (*emModel, error) {

	reportFit(phase, rounds, maxIters, tol)

	// 第一轮使用确定性的数据驱动种子（每序列背景相对得分最高的窗口），
	// 其余轮次为独立随机初始化；所有轮次共同竞争最大似然。
	// 未收敛的轮次直接丢弃；全部轮次都不收敛时返回 ErrEMNotConverged。
	var best *emModel
	if model, err := runEMSeeded(windows, width, bg, pseudocount, tol, maxIters,
		dataDrivenSeeds(windows, width, bg)); err == nil {
		best = model
	}
	failures := 0
	if best == nil {
		failures++
	}
	for round := 0; round < rounds; round++ {
		roundRand := streamRand(seed, baseStream+uint64(round))
		model, err := runEM(windows, width, bg, pseudocount, tol, maxIters, roundRand)
		if err != nil {
			failures++
			continue
		}
		if best == nil || model.loglik > best.loglik {
			best = model
		}
	}
	if best == nil {
		return nil, fmt.Errorf("%w: all %d init rounds (1 deterministic + %d random) hit the iteration limit",
			ErrEMNotConverged, failures, rounds)
	}
	return best, nil
}

// dataDrivenSeeds 为每条序列选出在背景下最“意外”的窗口作为确定性初始位点。
func dataDrivenSeeds(windows []window, width int, bg Background) []int {
	bySeq := make(map[int][]int)
	var seqOrder []int
	for wi, w := range windows {
		if _, ok := bySeq[w.seqIndex]; !ok {
			seqOrder = append(seqOrder, w.seqIndex)
		}
		bySeq[w.seqIndex] = append(bySeq[w.seqIndex], wi)
	}
	sort.Ints(seqOrder)

	seedIdx := make([]int, 0, len(seqOrder))
	for _, si := range seqOrder {
		bestWi := bySeq[si][0]
		bestScore := backgroundOdds(windows[bestWi].data, bg)
		for _, wi := range bySeq[si][1:] {
			score := backgroundOdds(windows[wi].data, bg)
			if score > bestScore {
				bestScore = score
				bestWi = wi
			}
		}
		seedIdx = append(seedIdx, bestWi)
	}
	return seedIdx
}

// backgroundOdds 返回窗口相对背景的总对数几率，作为无模体先验下的稀有度。
func backgroundOdds(window []byte, bg Background) float64 {
	score := 0.0
	for _, base := range window {
		p := bg.Prob(base)
		if p <= 0 {
			p = 0.25
		}
		score += math.Log(p / 0.25)
	}
	return score
}

// runEM 执行一次随机初始化的期望最大化。
func runEM(windows []window, width int, bg Background, pseudocount, tol float64,
	maxIters int, r *rng) (*emModel, error) {
	bySeq := make(map[int][]int)
	var seqOrder []int
	for wi, w := range windows {
		if _, ok := bySeq[w.seqIndex]; !ok {
			seqOrder = append(seqOrder, w.seqIndex)
		}
		bySeq[w.seqIndex] = append(bySeq[w.seqIndex], wi)
	}
	sort.Ints(seqOrder)

	initialIdx := make([]int, 0, len(seqOrder))
	for _, si := range seqOrder {
		candidates := bySeq[si]
		initialIdx = append(initialIdx, candidates[r.intn(len(candidates))])
	}
	return runEMSeeded(windows, width, bg, pseudocount, tol, maxIters, initialIdx)
}

// runEMSeeded 从给定的每序列初始窗口索引出发执行 EM。
// 迭代最多执行 maxIters 个 E/M 轮次；达到上限仍不满足
// |loglik-prevLoglik| <= tol 时返回包装 ErrEMNotConverged 的错误，
// 保证退化输入、近似重复输入与极端参数下也在有限步内结束，不返回半成品。
func runEMSeeded(windows []window, width int, bg Background, pseudocount, tol float64,
	maxIters int, initialIdx []int) (*emModel, error) {

	bySeq := make(map[int][]int)
	var seqOrder []int
	for wi, w := range windows {
		if _, ok := bySeq[w.seqIndex]; !ok {
			seqOrder = append(seqOrder, w.seqIndex)
		}
		bySeq[w.seqIndex] = append(bySeq[w.seqIndex], wi)
	}
	sort.Ints(seqOrder)

	initial := make([][]byte, 0, len(initialIdx))
	for _, wi := range initialIdx {
		initial = append(initial, windows[wi].data)
	}

	matrix, err := BuildPWM(initial, bg, pseudocount)
	if err != nil {
		// 初始片段来自合法窗口，构建不应失败。
		panic(err)
	}

	// 按序列组织窗口索引，保证 E 步在每条序列内归一化。
	weights := make([]float64, len(windows))
	prevLoglik := math.Inf(1)
	for iter := 0; iter < maxIters; iter++ {
		// E 步：在当前矩阵下计算每个候选窗口的后验权重（每序列内归一化）。
		loglik := 0.0
		for _, si := range seqOrder {
			candidates := bySeq[si]
			maxWeight := math.Inf(-1)
			for _, wi := range candidates {
				s := matrix.Score(windows[wi].data)
				weights[wi] = s
				if s > maxWeight {
					maxWeight = s
				}
			}
			sum := 0.0
			for _, wi := range candidates {
				weights[wi] = math.Exp(weights[wi] - maxWeight)
				sum += weights[wi]
			}
			for _, wi := range candidates {
				weights[wi] /= sum
			}
			loglik += maxWeight + math.Log(sum)
		}

		// NaN/Inf 防护：退化输入下似然可能被污染，此时立即以明确的不收敛
		// 结论结束，不继续迭代也不返回未定义的模型。
		if math.IsNaN(loglik) || math.IsInf(loglik, 0) {
			return nil, fmt.Errorf("%w: log-likelihood became non-finite at iteration %d",
				ErrEMNotConverged, iter+1)
		}
		if math.Abs(loglik-prevLoglik) <= tol {
			return selectSites(windows, weights, matrix, seqOrder, bySeq, loglik), nil
		}
		prevLoglik = loglik

		// M 步：按期望权重更新 PWM（加权计数 + 背景伪计数）。
		matrix = weightedPWM(windows, weights, width, bg, pseudocount)
	}
	return nil, fmt.Errorf("%w after %d iterations (tol=%g)",
		ErrEMNotConverged, maxIters, tol)
}

// weightedPWM 根据期望权重构建加权 PWM。
func weightedPWM(windows []window, weights []float64, width int, bg Background, pseudocount float64) *PWM {
	logOdds := make([][]float64, width)
	totalWeight := 0.0
	for _, w := range weights {
		totalWeight += w
	}
	for pos := 0; pos < width; pos++ {
		var counts [4]float64
		for wi, w := range windows {
			counts[baseIndex(w.data[pos])] += weights[wi]
		}
		row := make([]float64, 4)
		for b := 0; b < 4; b++ {
			prior := bg.Freq[b]
			if prior <= 0 {
				prior = 0.25
			}
			freq := (counts[b] + pseudocount*prior) / (totalWeight + pseudocount)
			row[b] = math.Log(freq / prior)
		}
		logOdds[pos] = row
	}

	m := &PWM{Width: width, LogOdds: logOdds, bg: bg}
	for _, row := range logOdds {
		lo, hi := row[0], row[0]
		for _, v := range row[1:] {
			if v < lo {
				lo = v
			}
			if v > hi {
				hi = v
			}
		}
		m.min += lo
		m.max += hi
	}
	return m
}

// selectSites 从最终权重中按每条序列取最大后验窗口，并记录信息含量。
func selectSites(windows []window, weights []float64, matrix *PWM,
	seqOrder []int, bySeq map[int][]int, loglik float64) *emModel {

	sites := make([]Site, 0, len(seqOrder))
	for _, si := range seqOrder {
		bestWi := bySeq[si][0]
		bestW := weights[bestWi]
		for _, wi := range bySeq[si][1:] {
			if weights[wi] > bestW {
				bestW = weights[wi]
				bestWi = wi
			}
		}
		w := windows[bestWi]
		sites = append(sites, Site{
			SeqIndex: si,
			Offset:   w.offset,
			Score:    matrix.NormalizedScore(w.data),
		})
	}
	sort.Slice(sites, func(i, j int) bool {
		if sites[i].SeqIndex != sites[j].SeqIndex {
			return sites[i].SeqIndex < sites[j].SeqIndex
		}
		return sites[i].Offset < sites[j].Offset
	})
	return &emModel{matrix: matrix, sites: sites, loglik: loglik}
}
