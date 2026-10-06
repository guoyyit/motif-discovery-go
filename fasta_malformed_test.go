package motif

import (
	"errors"
	"strings"
	"testing"
)

// mustParseError 断言输入被以 *ParseError 拒绝，且位置与预期一致，
// 并且不返回任何部分解析出的序列集合。
func mustParseError(t *testing.T, input string, wantLine, wantCol int) {
	t.Helper()
	seqs, err := ParseFASTA(strings.NewReader(input))
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("expected *ParseError, got %T: %v", err, err)
	}
	if pe.Line != wantLine || pe.Column != wantCol {
		t.Fatalf("position = %d:%d, want %d:%d (input %q)", pe.Line, pe.Column, wantLine, wantCol, input)
	}
	if seqs != nil {
		t.Fatalf("malformed input must not yield partial sequences, got %d", len(seqs))
	}
}

func TestFASTAHeaderWithoutSequence(t *testing.T) {
	// 头之后直接是下一个头：错误定位在空记录的头行。
	mustParseError(t, ">s1\n>s2\nACGT\n", 1, 1)
	// 头之后直到文件结束都没有序列数据。
	mustParseError(t, ">s1\nACGT\n>s2\n", 3, 1)
	// 头之后只有注释与空行，仍属无序列数据。
	mustParseError(t, ">s1\n; comment\n\n>s2\nACGT\n", 1, 1)
}

func TestFASTAEmptyNameAndEarlyData(t *testing.T) {
	mustParseError(t, ">\nACGT\n", 1, 1)
	mustParseError(t, ">   \nACGT\n", 1, 1)
	mustParseError(t, "ACGT\n>s1\n", 1, 1)
	mustParseError(t, "\n\nACGT\n>s1\nACGT\n", 3, 1)
}

func TestFASTAIllegalBasePosition(t *testing.T) {
	mustParseError(t, ">s1\nACGT\nACGX\n", 3, 4)
	mustParseError(t, ">s1\nAC GT X\n", 2, 7)
}

func TestFASTAOverlongLine(t *testing.T) {
	// 小限额：第 2 行有 6 个碱基，限额 5，越界位置为列 6。
	_, err := parseFASTALimits(strings.NewReader(">s1\nACGTAC\n"), 5, 100)
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("expected *ParseError, got %v", err)
	}
	if pe.Line != 2 || pe.Column != 6 {
		t.Fatalf("position = %d:%d, want 2:6", pe.Line, pe.Column)
	}

	// 默认限额：超过 1 MiB 的单行必须被拒绝且带位置。
	long := strings.Repeat("A", maxFASTALineLength+1)
	seqs, err := ParseFASTA(strings.NewReader(">s1\n" + long + "\n"))
	if !errors.As(err, &pe) {
		t.Fatalf("expected *ParseError for overlong line, got %v", err)
	}
	if pe.Line != 2 || pe.Column != maxFASTALineLength+1 {
		t.Fatalf("position = %d:%d, want 2:%d", pe.Line, pe.Column, maxFASTALineLength+1)
	}
	if seqs != nil {
		t.Fatal("overlong line must not yield partial sequences")
	}

	// 恰好等于限额的行是合法的。
	ok := strings.Repeat("A", maxFASTALineLength)
	if _, err := ParseFASTA(strings.NewReader(">s1\n" + ok + "\n")); err != nil {
		t.Fatalf("line at exactly the limit must parse: %v", err)
	}
}

func TestFASTAOversizedRecord(t *testing.T) {
	// 限额 5 个碱基：第 2 行写入第 6 个碱基时越界，位置为该行列 2。
	_, err := parseFASTALimits(strings.NewReader(">s1\nAAAAA\nAC\n"), 100, 5)
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("expected *ParseError, got %v", err)
	}
	if pe.Line != 3 || pe.Column != 1 {
		t.Fatalf("position = %d:%d, want 3:1", pe.Line, pe.Column)
	}

	// 记录超限时不返回部分序列集合。
	seqs, err := parseFASTALimits(strings.NewReader(">s1\nAAAAAA\n>s2\nACGT\n"), 100, 5)
	if !errors.As(err, &pe) {
		t.Fatalf("expected *ParseError, got %v", err)
	}
	if seqs != nil {
		t.Fatal("oversized record must not yield partial sequences")
	}
}

func TestFASTANoTrailingNewline(t *testing.T) {
	// 最后一行没有换行符也必须正常解析。
	seqs, err := ParseFASTA(strings.NewReader(">s1\nACGT"))
	if err != nil {
		t.Fatal(err)
	}
	if len(seqs) != 1 || string(seqs[0].Data) != "ACGT" {
		t.Fatalf("unexpected parse result: %+v", seqs)
	}
	// 无换行符的畸形行同样要定位。
	mustParseError(t, ">s1\nACGX", 2, 4)
}
