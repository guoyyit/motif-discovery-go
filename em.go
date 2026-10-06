package motif

import (
	"errors"
	"fmt"
	"math"
	"sort"
)

// ErrNotConverged 表示 EM 在达到 MaxIterations 迭代上限后仍未满足收敛容差。
// 该错误可通过 errors.Is 区分；出现时不返回任何半成品模型。
var ErrNotConverged = errors.New("motif: EM did not converge within MaxIterations")

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
// 每轮 EM 最多执行 maxIter 次迭代；未在限内收敛的轮次被丢弃，
// 所有轮次均不收敛时返回 ErrNotConverged，绝不返回半成品模型。
func fitOneMotif(windows []window, width int, bg Background, pseudocount, tol float64,
	rounds, maxIter int, seed uint64, baseStream uint64) (*emModel, error) {

	// 第一轮使用确定性的数据驱动种子（每序列背景相对得分最高的窗口），
	// 其余轮次为独立随机初始化；所有轮次共同竞争最大似然。
	var best *emModel
	if model, ok := runEMSeeded(windows, width, bg, pseudocount, tol, maxIter,
		dataDrivenSeeds(windows, width, bg)); ok {
		best = model
	}
	for round := 0; round < rounds; round++ {
		roundRand := streamRand(seed, baseStream+uint64(round))
		model, ok := runEM(windows, width, bg, pseudocount, tol, maxIter, roundRand)
		if !ok {
			continue
		}
		if best == nil || model.loglik > best.loglik {
			best = model
		}
	}
	if best == nil {
		return nil, fmt.Errorf("%w (rounds=%d, maxIterations=%d, tol=%g)",
			ErrNotConverged, rounds+1, maxIter, tol)
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
// 返回的布尔值表示是否在 maxIter 步内收敛。
func runEM(windows []window, width int, bg Background, pseudocount, tol float64,
	maxIter int, r *rng) (*emModel, bool) {

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
	return runEMSeeded(windows, width, bg, pseudocount, tol, maxIter, initialIdx)
}

// runEMSeeded 从给定的每序列初始窗口索引出发执行 EM。
// 最多执行 maxIter 次 E 步；达到上限仍未满足收敛容差时返回 ok=false，
// 此时不返回任何模型，保证调用方不会拿到半成品结果。
func runEMSeeded(windows []window, width int, bg Background, pseudocount, tol float64,
	maxIter int, initialIdx []int) (*emModel, bool) {

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
	for iter := 0; iter < maxIter; iter++ {
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

		if math.Abs(loglik-prevLoglik) <= tol {
			return selectSites(windows, weights, matrix, seqOrder, bySeq, loglik), true
		}
		prevLoglik = loglik

		// M 步：按期望权重更新 PWM（加权计数 + 背景伪计数）。
		matrix = weightedPWM(windows, weights, width, bg, pseudocount)
	}
	// 达到迭代上限仍未收敛：不返回模型，由调用方决定丢弃或报错。
	return nil, false
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
