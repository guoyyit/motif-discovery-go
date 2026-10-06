package motif

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
)

// checkpoint 是任务在磁盘上的可恢复状态。
// 每个工作单元（发现第 k 个模体、第 b 次重排）完成后立即落盘：
// 恢复时直接复用已完成单元，未完成单元因使用独立随机流而可逐位重演，
// 因此不会产生部分结果或重复计算。
type checkpoint struct {
	ID        string     `json:"id"`
	Status    string     `json:"status"`
	Error     string     `json:"error,omitempty"`
	Config    Config     `json:"config"`
	Sequences []Sequence `json:"sequences"`

	// DiscoverDone 已完成的发现单元数；已发现模体的可序列化表示。
	DiscoverDone     int          `json:"discover_done"`
	DiscoverComplete bool         `json:"discover_complete"`
	SavedModels      []savedModel `json:"saved_models"`

	// ShuffleDone 已完成的重排次数；NullScores[k] 为每个模体累积的零分布。
	ShuffleDone int         `json:"shuffle_done"`
	NullScores  [][]float64 `json:"null_scores"`

	// Result 在任务完成后写入。
	Result *savedResult `json:"result,omitempty"`
}

type savedModel struct {
	Width   int         `json:"width"`
	LogOdds [][]float64 `json:"log_odds"`
	Bg      [4]float64  `json:"bg"`
	Sites   []Site      `json:"sites"`
	Loglik  float64     `json:"loglik"`
}

type savedResult struct {
	Motifs []savedMotif `json:"motifs"`
}

type savedMotif struct {
	Index          int         `json:"index"`
	Width          int         `json:"width"`
	LogOdds        [][]float64 `json:"log_odds"`
	Bg             [4]float64  `json:"bg"`
	Sites          []Site      `json:"sites"`
	Information    float64     `json:"information"`
	Statistic      float64     `json:"statistic"`
	NullCount      int         `json:"null_count"`
	PValue         float64     `json:"p_value"`
	AdjustedPValue float64     `json:"adjusted_p_value"`
	Significant    bool        `json:"significant"`
}

func checkpointPath(dir, id string) string {
	return filepath.Join(dir, id+".json")
}

func toSavedModel(m *emModel) savedModel {
	return savedModel{
		Width:   m.matrix.Width,
		LogOdds: m.matrix.LogOdds,
		Bg:      m.matrix.bg.Freq,
		Sites:   m.sites,
		Loglik:  m.loglik,
	}
}

func fromSavedModel(sm savedModel) *emModel {
	bg := Background{Freq: sm.Bg}
	matrix := &PWM{Width: sm.Width, LogOdds: sm.LogOdds, bg: bg}
	for _, row := range sm.LogOdds {
		lo, hi := row[0], row[0]
		for _, v := range row[1:] {
			if v < lo {
				lo = v
			}
			if v > hi {
				hi = v
			}
		}
		matrix.min += lo
		matrix.max += hi
	}
	return &emModel{matrix: matrix, sites: sm.Sites, loglik: sm.Loglik}
}

func toSavedResult(r *Result) savedResult {
	out := savedResult{Motifs: make([]savedMotif, len(r.Motifs))}
	for i, m := range r.Motifs {
		out.Motifs[i] = savedMotif{
			Index:          m.Index,
			Width:          m.Matrix.Width,
			LogOdds:        m.Matrix.LogOdds,
			Bg:             m.Matrix.bg.Freq,
			Sites:          m.Sites,
			Information:    m.Information,
			Statistic:      m.Statistic,
			NullCount:      m.NullCount,
			PValue:         m.PValue,
			AdjustedPValue: m.AdjustedPValue,
			Significant:    m.Significant,
		}
	}
	return out
}

func fromSavedResult(sr savedResult) *Result {
	r := &Result{Motifs: make([]Motif, len(sr.Motifs))}
	for i, sm := range sr.Motifs {
		bg := Background{Freq: sm.Bg}
		matrix := &PWM{Width: sm.Width, LogOdds: sm.LogOdds, bg: bg}
		for _, row := range sm.LogOdds {
			lo, hi := row[0], row[0]
			for _, v := range row[1:] {
				if v < lo {
					lo = v
				}
				if v > hi {
					hi = v
				}
			}
			matrix.min += lo
			matrix.max += hi
		}
		r.Motifs[i] = Motif{
			Index:          sm.Index,
			Matrix:         matrix,
			Sites:          sm.Sites,
			Information:    sm.Information,
			Statistic:      sm.Statistic,
			NullCount:      sm.NullCount,
			PValue:         sm.PValue,
			AdjustedPValue: sm.AdjustedPValue,
			Significant:    sm.Significant,
		}
	}
	return r
}

func loadCheckpoint(dir, id string) (*checkpoint, error) {
	data, err := os.ReadFile(checkpointPath(dir, id))
	if err != nil {
		return nil, err
	}
	var cp checkpoint
	if err := json.Unmarshal(data, &cp); err != nil {
		return nil, fmt.Errorf("checkpoint %s: corrupt state: %w", id, err)
	}
	return &cp, nil
}

// store 原子写入检查点：先写临时文件再 rename，保证中断时状态完整。
func (cp *checkpoint) store(dir string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	data, err := json.MarshalIndent(cp, "", "  ")
	if err != nil {
		return err
	}
	final := checkpointPath(dir, cp.ID)
	tmp := final + ".tmp"
	if err := os.WriteFile(tmp, data, 0o644); err != nil {
		return err
	}
	return os.Rename(tmp, final)
}
