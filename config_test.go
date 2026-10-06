package motif

import (
	"math"
	"strings"
	"testing"
)

// TestValidateRejectsIllegalValues 显式给出的越界值必须直接报错，
// 说明字段名与合法范围，不得静默替换为默认值后继续运行。
func TestValidateRejectsIllegalValues(t *testing.T) {
	cases := map[string]struct {
		mutate func(*Config)
		field  string
	}{
		"negative motif length":   {func(c *Config) { c.MotifLength = -5 }, "MotifLength"},
		"negative expect motifs":  {func(c *Config) { c.ExpectMotifs = -1 }, "ExpectMotifs"},
		"negative init rounds":    {func(c *Config) { c.InitRounds = -2 }, "InitRounds"},
		"negative max iterations": {func(c *Config) { c.MaxIterations = -10 }, "MaxIterations"},
		"negative shuffles":       {func(c *Config) { c.Shuffles = -3 }, "Shuffles"},
		"negative convergence":    {func(c *Config) { c.ConvergenceTol = -1e-3 }, "ConvergenceTol"},
		"nan convergence":         {func(c *Config) { c.ConvergenceTol = math.NaN() }, "ConvergenceTol"},
		"negative score threshold": {
			func(c *Config) { c.ScoreThreshold = -0.5 }, "ScoreThreshold"},
		"score threshold above 1": {func(c *Config) { c.ScoreThreshold = 1.5 }, "ScoreThreshold"},
		"nan score threshold":     {func(c *Config) { c.ScoreThreshold = math.NaN() }, "ScoreThreshold"},
		"negative alpha":          {func(c *Config) { c.SignificanceAlpha = -0.1 }, "SignificanceAlpha"},
		"alpha above 1":           {func(c *Config) { c.SignificanceAlpha = 1.5 }, "SignificanceAlpha"},
		"nan alpha":               {func(c *Config) { c.SignificanceAlpha = math.NaN() }, "SignificanceAlpha"},
		"negative pseudocount":    {func(c *Config) { c.Pseudocount = -0.5 }, "Pseudocount"},
		"conflicting mask flags": {
			func(c *Config) { c.MaskLowComplexity, c.NoMaskLowComplexity = true, true }, "Mask"},
	}
	for name, tc := range cases {
		t.Run(name, func(t *testing.T) {
			cfg := testConfig(8)
			tc.mutate(&cfg)
			if _, err := cfg.Validate(); err == nil {
				t.Fatalf("expected error for %s", name)
			} else if !strings.Contains(err.Error(), tc.field) {
				t.Fatalf("error must name the field %q, got: %v", tc.field, err)
			}
			// Discover 入口同样必须拒绝，而不是替换默认值后继续。
			seqs := []Sequence{{Name: "s", Data: []byte("ACGTACGTACGTACGT")}}
			if _, err := Discover(seqs, cfg); err == nil {
				t.Fatalf("Discover accepted illegal config for %s", name)
			}
		})
	}
}

// TestValidateFillsDefaults 零值字段视为未配置，填入默认值。
func TestValidateFillsDefaults(t *testing.T) {
	cfg, err := (Config{}).Validate()
	if err != nil {
		t.Fatalf("zero config must be valid: %v", err)
	}
	if cfg.MotifLength != defaultMotifLength ||
		cfg.ExpectMotifs != defaultExpectMotifs ||
		cfg.InitRounds != defaultInitRounds ||
		cfg.ConvergenceTol != defaultConvergenceTol ||
		cfg.MaxIterations != defaultMaxIterations ||
		cfg.ScoreThreshold != defaultScoreThreshold ||
		cfg.SignificanceAlpha != defaultAlpha ||
		cfg.Shuffles != defaultShuffles ||
		cfg.Pseudocount != defaultPseudocount {
		t.Fatalf("defaults not filled: %+v", cfg)
	}
}

// TestValidateKeepsExplicitConfig 已有的显式合法配置必须按原语义生效，
// 校验不得改动任何显式字段。
func TestValidateKeepsExplicitConfig(t *testing.T) {
	explicit := Config{
		Seed:              99,
		MotifLength:       12,
		ExpectMotifs:      2,
		InitRounds:        7,
		ConvergenceTol:    1e-5,
		MaxIterations:     123,
		ScoreThreshold:    0.9,
		SignificanceAlpha: 0.01,
		Shuffles:          77,
		Pseudocount:       0.25,
		MaskLowComplexity: true,
	}
	got, err := explicit.Validate()
	if err != nil {
		t.Fatalf("explicit config rejected: %v", err)
	}
	if got != explicit {
		t.Fatalf("explicit config altered:\nwant %+v\ngot  %+v", explicit, got)
	}

	// 显式关闭掩码的旧行为逃生口同样必须保留。
	open := Config{MotifLength: 8, NoMaskLowComplexity: true}
	got, err = open.Validate()
	if err != nil {
		t.Fatalf("NoMaskLowComplexity config rejected: %v", err)
	}
	if !got.NoMaskLowComplexity || maskEnabled(got) {
		t.Fatalf("explicit NoMaskLowComplexity not honored: %+v", got)
	}
}

// TestMaskDefaultOn 默认（零值）配置必须开启低复杂度掩码，
// 与显式 MaskLowComplexity=true 的结论一致。
func TestMaskDefaultOn(t *testing.T) {
	def, err := (Config{}).Validate()
	if err != nil {
		t.Fatal(err)
	}
	if !maskEnabled(def) {
		t.Fatal("masking must be on by default")
	}
	legacy, err := (Config{MaskLowComplexity: true}).Validate()
	if err != nil {
		t.Fatal(err)
	}
	if !maskEnabled(legacy) {
		t.Fatal("explicit MaskLowComplexity=true must keep masking on")
	}
}
