package motif

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strings"
)

// Sequence 表示一条带名称的核酸序列。
type Sequence struct {
	Name        string
	Description string
	Data        []byte
}

// ParseError 记录 FASTA 解析错误的位置（从 1 开始计数）。
type ParseError struct {
	Line   int
	Column int
	Msg    string
}

func (e *ParseError) Error() string { return e.Msg }

// ParseFASTA 解析 FASTA 格式的序列集合。
//
// 规则：
//   - 以 '>' 开头的行开启一条新记录，其后第一个空白分隔字段为名称，其余为描述。
//   - 以 ';' 或 '#' 开头的行以及空行视为注释，可出现在任意位置。
//   - 序列允许大小写混用，解析后统一归一化为大写；允许的碱基字符为 A/C/G/T/N，
//     其中 N 表示未知碱基（低复杂度掩码同样使用 N）。序列行中的空格与制表符会被忽略。
//   - 任何非法字符或格式问题都会返回 *ParseError，注明从 1 开始的行号与列号。
func ParseFASTA(r io.Reader) ([]Sequence, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)

	var sequences []Sequence
	var current *Sequence
	var data bytes.Buffer

	flush := func() {
		if current != nil {
			current.Data = append(current.Data, data.Bytes()...)
			sequences = append(sequences, *current)
			data.Reset()
			current = nil
		}
	}

	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := scanner.Text()
		trimmed := strings.TrimSpace(line)

		switch {
		case trimmed == "":
			continue
		case trimmed[0] == ';' || trimmed[0] == '#':
			continue
		case trimmed[0] == '>':
			flush()
			header := strings.TrimSpace(trimmed[1:])
			if header == "" {
				return nil, &ParseError{Line: lineNo, Column: 1,
					Msg: fmt.Sprintf("fasta: line %d column 1: header line '>' is missing a sequence name", lineNo)}
			}
			fields := strings.Fields(header)
			current = &Sequence{
				Name:        fields[0],
				Description: strings.TrimSpace(strings.TrimPrefix(header, fields[0])),
			}
		default:
			if current == nil {
				return nil, &ParseError{Line: lineNo, Column: 1,
					Msg: fmt.Sprintf("fasta: line %d column 1: sequence data appears before any '>' header", lineNo)}
			}
			for col, ch := range []byte(line) {
				if ch == ' ' || ch == '\t' {
					continue
				}
				upper := ch
				if upper >= 'a' && upper <= 'z' {
					upper -= 'a' - 'A'
				}
				switch upper {
				case 'A', 'C', 'G', 'T', 'N':
					data.WriteByte(upper)
				default:
					return nil, &ParseError{Line: lineNo, Column: col + 1,
						Msg: fmt.Sprintf("fasta: line %d column %d: illegal character %q in sequence %q",
							lineNo, col+1, rune(ch), current.Name)}
				}
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("fasta: read error: %w", err)
	}
	flush()

	if len(sequences) == 0 {
		return nil, &ParseError{Line: 0, Column: 0, Msg: "fasta: input contains no sequences"}
	}
	return sequences, nil
}
