package motif

import (
	"errors"
	"strings"
	"testing"
)

func TestParseFASTA(t *testing.T) {
	input := strings.NewReader("; comment line\n" +
		"# another comment\n" +
		">seq1 first sequence\n" +
		"acgt\n" +
		"ACGT Nn\n" +
		"\n" +
		">seq2\n" +
		"tttt\n")

	seqs, err := ParseFASTA(input)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(seqs) != 2 {
		t.Fatalf("got %d sequences, want 2", len(seqs))
	}
	if seqs[0].Name != "seq1" || seqs[0].Description != "first sequence" {
		t.Fatalf("unexpected header: %+v", seqs[0])
	}
	if string(seqs[0].Data) != "ACGTACGTNN" {
		t.Fatalf("unexpected data: %q", seqs[0].Data)
	}
	if string(seqs[1].Data) != "TTTT" {
		t.Fatalf("unexpected seq2 data: %q", seqs[1].Data)
	}
}

func TestParseFASTAErrors(t *testing.T) {
	cases := map[string]string{
		"illegal char": ">s1\nacgx\n",
		"data first":   "acgt\n>s1\n",
		"empty header": ">\n",
		"no sequences": "; only comment\n",
	}
	for name, in := range cases {
		t.Run(name, func(t *testing.T) {
			_, err := ParseFASTA(strings.NewReader(in))
			var pe *ParseError
			if !errors.As(err, &pe) {
				t.Fatalf("expected *ParseError, got %v", err)
			}
		})
	}

	_, err := ParseFASTA(strings.NewReader(">s1\nacgx\n"))
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("expected *ParseError, got %v", err)
	}
	if pe.Line != 2 || pe.Column != 4 {
		t.Fatalf("unexpected position line=%d col=%d", pe.Line, pe.Column)
	}
}
