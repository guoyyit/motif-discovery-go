package motif

import (
	"strings"
	"testing"
)

func TestFASTAToDiscover(t *testing.T) {
	var sb strings.Builder
	motif := "ACGTATGC"
	r := streamRand(424, globalStream(9002, 0))
	bases := []byte("ACGT")
	for i := 0; i < 15; i++ {
		sb.WriteString(">seq")
		sb.WriteByte(byte('a' + i%26))
		sb.WriteString(" sample record\n")
		seq := make([]byte, 50)
		for j := range seq {
			seq[j] = bases[r.intn(4)]
		}
		pos := r.intn(50 - len(motif) + 1)
		copy(seq[pos:pos+len(motif)], motif)
		sb.Write(seq)
		sb.WriteByte('\n')
	}

	seqs, err := ParseFASTA(strings.NewReader(sb.String()))
	if err != nil {
		t.Fatal(err)
	}
	cfg := testConfig(len(motif))
	cfg.Shuffles = 25
	result, err := Discover(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(result.Motifs) != 1 || !result.Motifs[0].Significant {
		t.Fatalf("expected one significant motif, got %+v", result.Motifs)
	}
	if result.Motifs[0].Matrix.Consensus() != motif {
		t.Fatalf("consensus=%q, want %q", result.Motifs[0].Matrix.Consensus(), motif)
	}
}

func TestParseErrorPosition(t *testing.T) {
	_, err := ParseFASTA(strings.NewReader(">s\nACGT\n1234\n"))
	pe, ok := err.(*ParseError)
	if !ok {
		t.Fatalf("expected *ParseError, got %T", err)
	}
	if pe.Line != 3 || pe.Column != 1 {
		t.Fatalf("error position line=%d col=%d, want 3:1", pe.Line, pe.Column)
	}
	if !strings.Contains(pe.Error(), "line 3") {
		t.Fatalf("error message should mention line: %q", pe.Error())
	}
}
