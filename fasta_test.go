package motif

import "testing"

func TestParseFASTABasic(t *testing.T) {
	in := "; a comment line\n" +
		">seq1 first sequence\n" +
		"acgt\n" +
		"ACGT-tt\n" +
		"\n" +
		"  ; indented comment\n" +
		">seq2\n" +
		"ggg ccc\n"
	_, err := ParseFASTAStr(in)
	if err == nil {
		t.Fatalf("expected error for '-' character, got nil")
	}
	pe, ok := err.(*ParseError)
	if !ok {
		t.Fatalf("expected *ParseError, got %T", err)
	}
	if pe.Line != 4 || pe.Column != 5 {
		t.Fatalf("unexpected position line=%d col=%d", pe.Line, pe.Column)
	}
}

func TestParseFASTAOK(t *testing.T) {
	in := ">alpha gene one\nacgt\nTTTT\n> beta spaced name\ng g\n"
	seqs, err := ParseFASTAStr(in)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(seqs) != 2 {
		t.Fatalf("want 2 seqs, got %d", len(seqs))
	}
	if seqs[0].Name != "alpha" || seqs[0].Desc != "gene one" || string(seqs[0].Data) != "ACGTTTTT" {
		t.Fatalf("seq0 mismatch: %+v data=%q", seqs[0], seqs[0].Data)
	}
	if seqs[1].Name != "beta" || seqs[1].Desc != "spaced name" || string(seqs[1].Data) != "GG" {
		t.Fatalf("seq1 mismatch: %+v data=%q", seqs[1], seqs[1].Data)
	}
}

func TestParseFASTAErrors(t *testing.T) {
	cases := map[string]string{
		"data before header": "acgt\n>s\n",
		"empty record":       ">s\n>s2\nacgt\n",
		"empty input":        "; only comment\n",
		"illegal N":          ">s\nacgnt\n",
		"missing name":       ">\nacgt\n",
	}
	for name, in := range cases {
		if _, err := ParseFASTAStr(in); err == nil {
			t.Errorf("%s: expected error", name)
		}
	}
}
