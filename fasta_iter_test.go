package motif

import (
	"errors"
	"strings"
	"testing"
)

func parseErrAt(t *testing.T, input string, wantLine, wantCol int) {
	t.Helper()
	seqs, err := ParseFASTA(strings.NewReader(input))
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("want *ParseError at %d:%d, got %v (seqs=%d)", wantLine, wantCol, err, len(seqs))
	}
	if pe.Line != wantLine || pe.Column != wantCol {
		t.Fatalf("error at %d:%d, want %d:%d (%s)", pe.Line, pe.Column, wantLine, wantCol, pe.Msg)
	}
	if seqs != nil {
		t.Fatalf("malformed input must not return a partial sequence set, got %d", len(seqs))
	}
}

func TestFASTAEmptyHeaderPosition(t *testing.T) {
	// 空名称头（含只有空白的头）。
	parseErrAt(t, ">ok\nACGT\n>\nTT", 3, 1)
	parseErrAt(t, ">   \nACGT\n", 1, 1)
}

func TestFASTAHeaderWithoutSequence(t *testing.T) {
	// 文件结束时只有头没有序列数据：定位到该头所在行。
	parseErrAt(t, ">s1\nACGT\n>s2\n", 3, 1)
	// 下一条头出现前没有任何序列数据：定位到空头所在行（中间隔注释/空行）。
	parseErrAt(t, ">s1\n# comment\n\n>s2\nACGT\n", 1, 1)
}

func TestFASTADataBeforeHeader(t *testing.T) {
	parseErrAt(t, "ACGT\n>s1\n", 1, 1)
	parseErrAt(t, "; c\n\nTTTT\n>s\n", 3, 1)
}

func TestFASTAILlegalBasePosition(t *testing.T) {
	parseErrAt(t, ">s\nACGT\nACG X\n", 3, 5)
	// 列号按物理字节计数：'A'=1,'C'=2,'\t'=3,'G'=4,' '=5,'Z'=6。
	parseErrAt(t, ">s\nAC\tG Z\n", 2, 6)
}

func TestFASTAWhitespaceColumnsIgnored(t *testing.T) {
	seqs, err := ParseFASTA(strings.NewReader(">s\nAC GT\nTT\tNN\n"))
	if err != nil {
		t.Fatal(err)
	}
	if string(seqs[0].Data) != "ACGTTTNN" {
		t.Fatalf("unexpected data %q", seqs[0].Data)
	}
}

func TestFASTAOverLongLine(t *testing.T) {
	var sb strings.Builder
	sb.WriteString(">s\n")
	sb.WriteString(strings.Repeat("A", MaxFASTALineBytes+1))
	sb.WriteString("\nACGT\n")
	parseErrAt(t, sb.String(), 2, MaxFASTALineBytes+1)

	// 超长行之后的合法行不得造成行号错位（这里错误在超长行即终止）。
	seqs, err := ParseFASTA(strings.NewReader(">s\n" +
		strings.Repeat("A", MaxFASTALineBytes+1) + "\n>s2\nACGT\n"))
	if err == nil {
		t.Fatalf("expected over-long line error, got seqs=%d", len(seqs))
	}
}

func TestFASTAOverLargeRecord(t *testing.T) {
	// 多条未超长的行累积成超大记录：必须在越界字节处定位报错。
	var sb strings.Builder
	sb.WriteString(">big\n")
	// 每行 1000 个碱基，写到刚好越界为止。
	line := strings.Repeat("C", 1000)
	for written := 0; written <= MaxFASTARecordBytes; written += 1000 {
		sb.WriteString(line)
		sb.WriteByte('\n')
	}
	seqs, err := ParseFASTA(strings.NewReader(sb.String()))
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("want *ParseError for oversized record, got %v (seqs=%d)", err, len(seqs))
	}
	if pe.Line < 1 || pe.Column < 1 {
		t.Fatalf("error must carry 1-based line/column, got %d:%d", pe.Line, pe.Column)
	}
	if seqs != nil {
		t.Fatal("oversized record must not return partial sequences")
	}
}

func TestFASTANoPartialSetOnLaterBadRecord(t *testing.T) {
	// 第一条记录合法、第二条含非法碱基：整体失败，不得只返回第一条。
	parseErrAt(t, ">s1\nACGT\n>s2\nACGX\n", 4, 4)
}
