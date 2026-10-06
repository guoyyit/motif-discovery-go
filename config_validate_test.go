package motif

import (
	"math"
	"strings"
	"testing"
)

// TestValidateDefaultsFilled 缺省字段必须填入既定默认值。
func TestValidateDefaultsFilled(t *testing.T) {
	vc, err := Config{}.Validate()
	if err != nil {
		t.Fatalf("zero config should validate with defaults: %v", err)
	}
	if vc.MotifLength != defaultMotifLength ||
		vc.ExpectMotifs != defaultExpectMotifs ||
		vc.InitRounds != defaultInitRounds ||
		vc.ConvergenceTol != defaultConvergenceTol ||
		*vc.ScoreThreshold != defaultScoreThreshold ||
		*vc.SignificanceAlpha != defaultAlpha ||
		vc.Shuffles != defaultShuffles ||
		vc.Pseudocount != defaultPseudocount ||
		vc.MaxEMIterations != defaultMaxEMIterations {
		t.Fatalf("unexpected defaults: %+v", vc)
	}
	if !vc.maskLowComplexity() {
		t.Fatal("low-complexity masking must default to enabled")
	}
}

// TestValidateExplicitIllegalRejected 显式给出的越界值必须直接报错，
// 错误信息需包含字段名与取值范围，而不是被静默替换为默认值。
func TestValidateExplicitIllegalRejected(t *testing.T) {
	bad := []struct {
		name   string
		mutate func(*Config)
		field  string
	}{
		{"score threshold zero", func(c *Config) { c.ScoreThreshold = Float64(0) }, "ScoreThreshold"},
		{"score threshold negative", func(c *Config) { c.ScoreThreshold = Float64(-0.1) }, "ScoreThreshold"},
		{"score threshold > 1", func(c *Config) { c.ScoreThreshold = Float64(1.5) }, "ScoreThreshold"},
		{"score threshold NaN", func(c *Config) { c.ScoreThreshold = Float64(math.NaN()) }, "ScoreThreshold"},
		{"alpha negative", func(c *Config) { c.SignificanceAlpha = Float64(-0.05) }, "SignificanceAlpha"},
		{"alpha > 1", func(c *Config) { c.SignificanceAlpha = Float64(1.2) }, "SignificanceAlpha"},
		{"alpha zero", func(c *Config) { c.SignificanceAlpha = Float64(0) }, "SignificanceAlpha"},
		{"alpha NaN", func(c *Config) { c.SignificanceAlpha = Float64(math.NaN()) }, "SignificanceAlpha"},
		{"negative init rounds", func(c *Config) { c.InitRounds = -1 }, "InitRounds"},
		{"negative shuffles", func(c *Config) { c.Shuffles = -10 }, "Shuffles"},
		{"negative expect motifs", func(c *Config) { c.ExpectMotifs = -2 }, "ExpectMotifs"},
		{"negative motif length", func(c *Config) { c.MotifLength = -8 }, "MotifLength"},
		{"motif length 3", func(c *Config) { c.MotifLength = 3 }, "MotifLength"},
		{"negative max iter", func(c *Config) { c.MaxEMIterations = -1 }, "MaxEMIterations"},
		{"negative tol", func(c *Config) { c.ConvergenceTol = -1 }, "ConvergenceTol"},
		{"tol NaN", func(c *Config) { c.ConvergenceTol = math.NaN() }, "ConvergenceTol"},
		{"negative pseudocount", func(c *Config) { c.Pseudocount = -0.5 }, "Pseudocount"},
		{"pseudocount NaN", func(c *Config) { c.Pseudocount = math.NaN() }, "Pseudocount"},
	}
	for _, tc := range bad {
		t.Run(tc.name, func(t *testing.T) {
			cfg := testConfig(8)
			tc.mutate(&cfg)
			vc, err := cfg.Validate()
			if err == nil {
				t.Fatalf("expected validation error, got validated config: %+v", vc)
			}
			if !strings.Contains(err.Error(), tc.field) {
				t.Fatalf("error %q must name field %q", err.Error(), tc.field)
			}
			// 不得静默替换为默认值后继续运行：非法字段必须原样保留在报错信息中。
		})
	}
}

// TestValidateExplicitLegalPreserved 显式合法配置必须原样保留（原语义不被破坏）。
func TestValidateExplicitLegalPreserved(t *testing.T) {
	cfg := testConfig(8)
	cfg.MaxEMIterations = 7
	cfg.MaskLowComplexity = Bool(false)
	cfg.ScoreThreshold = Float64(0.42)
	cfg.SignificanceAlpha = Float64(0.2)
	cfg.ConvergenceTol = 1e-7
	cfg.Pseudocount = 1.25
	vc, err := cfg.Validate()
	if err != nil {
		t.Fatal(err)
	}
	if vc.MaxEMIterations != 7 || vc.maskLowComplexity() ||
		*vc.ScoreThreshold != 0.42 || *vc.SignificanceAlpha != 0.2 ||
		vc.ConvergenceTol != 1e-7 || vc.Pseudocount != 1.25 {
		t.Fatalf("explicit legal values were altered: %+v", vc)
	}
}

// TestValidateMaskingTristate 缺省与显式 true 等价，显式 false 关闭保护。
func TestValidateMaskingTristate(t *testing.T) {
	def, err := Config{}.Validate()
	if err != nil {
		t.Fatal(err)
	}
	on, err := Config{MaskLowComplexity: Bool(true)}.Validate()
	if err != nil {
		t.Fatal(err)
	}
	off, err := Config{MaskLowComplexity: Bool(false)}.Validate()
	if err != nil {
		t.Fatal(err)
	}
	if !def.maskLowComplexity() || !on.maskLowComplexity() || off.maskLowComplexity() {
		t.Fatalf("mask tristate wrong: default=%v explicit-on=%v explicit-off=%v",
			def.maskLowComplexity(), on.maskLowComplexity(), off.maskLowComplexity())
	}
}
