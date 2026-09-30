package motif

// shuffleBases 对序列中的已知碱基（A/C/G/T）执行确定性的 Fisher-Yates 重排。
// 未知碱基 N 保持在原位，因此重排后低复杂度掩码结构不变，
// 只打乱碱基排列以构成零假设（“无模体结构，仅有相同碱基组成”）。
// 重排在输入切片上原地进行，调用方需要保留原序列时应先拷贝。
func shuffleBases(seq []byte, r *rng) {
	var knownIdx []int
	for i, base := range seq {
		if baseIndex(base) >= 0 {
			knownIdx = append(knownIdx, i)
		}
	}
	for i := len(knownIdx) - 1; i > 0; i-- {
		j := r.intn(i + 1)
		a, b := knownIdx[i], knownIdx[j]
		seq[a], seq[b] = seq[b], seq[a]
		knownIdx[i], knownIdx[j] = knownIdx[j], knownIdx[i]
	}
}

// shuffleAll 返回所有序列的重排副本，保留输入不变。
func shuffleAll(seqs [][]byte, r *rng) [][]byte {
	out := make([][]byte, len(seqs))
	for i, seq := range seqs {
		dup := append([]byte(nil), seq...)
		shuffleBases(dup, r)
		out[i] = dup
	}
	return out
}
