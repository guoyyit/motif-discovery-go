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

// FASTA 解析的容量上限：任一限制被突破都会以带位置的 ParseError 拒绝，
// 避免畸形或恶意输入造成无界内存消耗。
const (
	// maxFASTALineLength 是单行允许的最大字节数（不含换行符）。
	maxFASTALineLength = 1 << 20 // 1 MiB
	// maxFASTARecordBases 是单条记录允许的最大碱基数。
	maxFASTARecordBases = 1 << 26 // 64 Mi bases
)

// ParseFASTA 解析 FASTA 格式的序列集合。
//
// 规则：
//   - 以 '>' 开头的行开启一条新记录，其后第一个空白分隔字段为名称，其余为描述。
//   - 以 ';' 或 '#' 开头的行以及空行视为注释，可出现在任意位置。
//   - 序列允许大小写混用，解析后统一归一化为大写；允许的碱基字符为 A/C/G/T/N，
//     其中 N 表示未知碱基（低复杂度掩码同样使用 N）。序列行中的空格与制表符会被忽略。
//   - 任何非法字符或格式问题都会返回 *ParseError，注明从 1 开始的行号与列号；
//     出错时不返回任何部分解析出的序列。
//
// 被拒绝的畸形输入包括：空名称头、只有头没有序列数据、序列数据出现在任何
// 头之前、非法碱基字符、超过 maxFASTALineLength 的超长单行，以及超过
// maxFASTARecordBases 的超大记录。
func ParseFASTA(r io.Reader) ([]Sequence, error) {
	return parseFASTALimits(r, maxFASTALineLength, maxFASTARecordBases)
}

// parseFASTALimits 是带显式容量上限的解析实现，便于测试小限额场景。
func parseFASTALimits(r io.Reader, maxLine, maxRecord int) ([]Sequence, error) {
	br := bufio.NewReader(r)

	var sequences []Sequence
	var current *Sequence
	headerLine := 0
	var data bytes.Buffer

	// flush 结束当前记录；只有头没有任何序列数据属于畸形输入。
	flush := func() error {
		if current == nil {
			return nil
		}
		if data.Len() == 0 {
			return &ParseError{Line: headerLine, Column: 1,
				Msg: fmt.Sprintf("fasta: line %d column 1: record %q has a header but no sequence data",
					headerLine, current.Name)}
		}
		current.Data = append(current.Data, data.Bytes()...)
		sequences = append(sequences, *current)
		data.Reset()
		current = nil
		return nil
	}

	lineNo := 0
	for {
		raw, readErr := br.ReadBytes('\n')
		if len(raw) > 0 {
			lineNo++
			line := raw
			if line[len(line)-1] == '\n' {
				line = line[:len(line)-1]
			}
			if len(line) > maxLine {
				return nil, &ParseError{Line: lineNo, Column: maxLine + 1,
					Msg: fmt.Sprintf("fasta: line %d column %d: line exceeds maximum length of %d bytes",
						lineNo, maxLine+1, maxLine)}
			}
			if err := parseLine(line, lineNo, &current, &headerLine, &data, maxRecord, flush); err != nil {
				return nil, err
			}
		}
		if readErr != nil {
			if readErr == io.EOF {
				break
			}
			return nil, fmt.Errorf("fasta: read error: %w", readErr)
		}
	}
	if err := flush(); err != nil {
		return nil, err
	}

	if len(sequences) == 0 {
		return nil, &ParseError{Line: 0, Column: 0, Msg: "fasta: input contains no sequences"}
	}
	return sequences, nil
}

// parseLine 处理单行内容；列号按原始行字节位置从 1 开始计数。
func parseLine(line []byte, lineNo int, current **Sequence, headerLine *int,
	data *bytes.Buffer, maxRecord int, flush func() error) error {

	trimmed := strings.TrimSpace(string(line))
	switch {
	case trimmed == "":
		return nil
	case trimmed[0] == ';' || trimmed[0] == '#':
		return nil
	case trimmed[0] == '>':
		if err := flush(); err != nil {
			return err
		}
		header := strings.TrimSpace(trimmed[1:])
		if header == "" {
			return &ParseError{Line: lineNo, Column: 1,
				Msg: fmt.Sprintf("fasta: line %d column 1: header line '>' is missing a sequence name", lineNo)}
		}
		fields := strings.Fields(header)
		*current = &Sequence{
			Name:        fields[0],
			Description: strings.TrimSpace(strings.TrimPrefix(header, fields[0])),
		}
		*headerLine = lineNo
		return nil
	default:
		if *current == nil {
			return &ParseError{Line: lineNo, Column: 1,
				Msg: fmt.Sprintf("fasta: line %d column 1: sequence data appears before any '>' header", lineNo)}
		}
		for col, ch := range line {
			if ch == ' ' || ch == '\t' {
				continue
			}
			upper := ch
			if upper >= 'a' && upper <= 'z' {
				upper -= 'a' - 'A'
			}
			switch upper {
			case 'A', 'C', 'G', 'T', 'N':
				if data.Len() >= maxRecord {
					return &ParseError{Line: lineNo, Column: col + 1,
						Msg: fmt.Sprintf("fasta: line %d column %d: record %q exceeds maximum of %d bases",
							lineNo, col+1, (*current).Name, maxRecord)}
				}
				data.WriteByte(upper)
			default:
				return &ParseError{Line: lineNo, Column: col + 1,
					Msg: fmt.Sprintf("fasta: line %d column %d: illegal character %q in sequence %q",
						lineNo, col+1, rune(ch), (*current).Name)}
			}
		}
		return nil
	}
}
