package motif

import (
	"bytes"
	"fmt"
	"io"
	"strings"
)

// 容量边界（默认值）。超过任一上限的输入在越界字节所在的
// 行/列位置返回 *ParseError；不会继续读取到内存爆炸。
const (
	// MaxFASTALineBytes 是单个物理行允许的最大字节数（含换行符前全部内容）。
	MaxFASTALineBytes = 1 << 20 // 1 MiB
	// MaxFASTARecordBytes 是单条记录序列数据允许的最大字节数（忽略空白后）。
	MaxFASTARecordBytes = 64 << 20 // 64 MiB
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

// cappedLineReader 以行（'\n' 分隔）为单位读取，并对单行长度设硬上限。
// 超长时在越界字节处返回 errLineTooLong，且仍会消费完该行剩余字节，
// 保证后续行号统计准确且不会无限缓冲。
type cappedLineReader struct {
	r     io.Reader
	limit int

	tail []byte // 尚未遇到换行符的读缓冲
}

var errLineTooLong = fmt.Errorf("line too long")

// readLine 返回下一行（不含行尾）。第二个返回值表示是否消耗了一个行号：
// 读到换行分隔的行或无换行结尾的末行时均为 true；输入彻底结束时为 false。
func (c *cappedLineReader) readLine() (line []byte, consumedLine bool, err error) {
	// 若缓冲中已有完整行（上一轮多读的内容），先返回它。
	if i := bytes.IndexByte(c.tail, '\n'); i >= 0 {
		line = c.tail[:i]
		c.tail = c.tail[i+1:]
		return line, true, nil
	}
	for {
		tmp := make([]byte, 64*1024)
		n, readErr := c.r.Read(tmp)
		c.tail = append(c.tail, tmp[:n]...)

		// 追加后立即判超长：无论本 chunk 是否同时带换行，超限即报错，
		// 且 discardRestOfLine 只保留换行之后的内容，行号不丢。
		if len(c.tail) > c.limit {
			c.discardRestOfLine()
			return nil, true, errLineTooLong
		}
		if i := bytes.IndexByte(c.tail, '\n'); i >= 0 {
			line = c.tail[:i]
			c.tail = c.tail[i+1:]
			return line, true, nil
		}
		if readErr != nil {
			if readErr == io.EOF {
				if len(c.tail) == 0 {
					return nil, false, nil
				}
				line = c.tail
				c.tail = nil
				return line, true, nil
			}
			return nil, false, readErr
		}
	}
}

// discardRestOfLine 消费当前行剩余字节（不保留），保持行计数准确且不无限缓冲。
func (c *cappedLineReader) discardRestOfLine() {
	for {
		if i := bytes.IndexByte(c.tail, '\n'); i >= 0 {
			c.tail = c.tail[i+1:]
			return
		}
		c.tail = c.tail[:0]
		tmp := make([]byte, 64*1024)
		n, err := c.r.Read(tmp)
		c.tail = append(c.tail, tmp[:n]...)
		if err != nil {
			c.tail = nil
			return
		}
	}
}

// ParseFASTA 解析 FASTA 格式的序列集合。
//
// 规则：
//   - 以 '>' 开头的行开启一条新记录，其后第一个空白分隔字段为名称，其余为描述。
//   - 以 ';' 或 '#' 开头的行以及空行视为注释，可出现在任意位置。
//   - 序列允许大小写混用，解析后统一归一化为大写；允许的碱基字符为 A/C/G/T/N，
//     其中 N 表示未知碱基（低复杂度掩码同样使用 N）。序列行中的空格与制表符会被忽略。
//   - 任何非法字符或格式问题（空名称头、只有头没有序列、数据出现在头之前、
//     非法碱基、超长单行、超大记录）都会返回 *ParseError，注明从 1 开始的
//     行号与列号。任一记录非法时整体失败，绝不返回部分序列集合。
func ParseFASTA(r io.Reader) ([]Sequence, error) {
	cl := &cappedLineReader{r: r, limit: MaxFASTALineBytes}

	var sequences []Sequence
	var current *Sequence
	currentLine := 0
	var data bytes.Buffer
	recordBytes := 0

	flush := func() bool {
		if current != nil {
			if recordBytes == 0 {
				return false
			}
			current.Data = append(current.Data, data.Bytes()...)
			sequences = append(sequences, *current)
			data.Reset()
			recordBytes = 0
			current = nil
		}
		return true
	}

	lineNo := 0
	for {
		raw, hasLine, err := cl.readLine()
		if err != nil {
			if err == errLineTooLong {
				// 超长行可能是文件最后一行（无换行），仍占用一个行号。
				lineNo++
				col := MaxFASTALineBytes + 1
				return nil, &ParseError{Line: lineNo, Column: col,
					Msg: fmt.Sprintf("fasta: line %d column %d: line exceeds maximum length of %d bytes",
						lineNo, col, MaxFASTALineBytes)}
			}
			return nil, fmt.Errorf("fasta: read error: %w", err)
		}
		if !hasLine {
			break
		}
		lineNo++
		line := raw
		if len(line) > 0 && line[len(line)-1] == '\r' {
			line = line[:len(line)-1]
		}
		trimmed := strings.TrimSpace(string(line))

		switch {
		case trimmed == "":
			continue
		case trimmed[0] == ';' || trimmed[0] == '#':
			continue
		case trimmed[0] == '>':
			if !flush() {
				return nil, &ParseError{Line: currentLine, Column: 1,
					Msg: fmt.Sprintf("fasta: line %d column 1: header %q has no sequence data before the next record",
						currentLine, current.Name)}
			}
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
			currentLine = lineNo
		default:
			if current == nil {
				return nil, &ParseError{Line: lineNo, Column: 1,
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
					if recordBytes >= MaxFASTARecordBytes {
						return nil, &ParseError{Line: lineNo, Column: col + 1,
							Msg: fmt.Sprintf("fasta: line %d column %d: record %q exceeds maximum size of %d bytes",
								lineNo, col+1, current.Name, MaxFASTARecordBytes)}
					}
					recordBytes++
					data.WriteByte(upper)
				default:
					return nil, &ParseError{Line: lineNo, Column: col + 1,
						Msg: fmt.Sprintf("fasta: line %d column %d: illegal character %q in sequence %q",
							lineNo, col+1, rune(ch), current.Name)}
				}
			}
		}
	}

	// 文件结束时最后一条记录仍打开：必须含至少一个碱基。
	if current != nil {
		if recordBytes == 0 {
			return nil, &ParseError{Line: currentLine, Column: 1,
				Msg: fmt.Sprintf("fasta: line %d column 1: header %q has no sequence data",
					currentLine, current.Name)}
		}
		flush()
	}

	if len(sequences) == 0 {
		return nil, &ParseError{Line: 0, Column: 0, Msg: "fasta: input contains no sequences"}
	}
	return sequences, nil
}
