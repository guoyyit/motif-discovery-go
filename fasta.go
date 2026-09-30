package motif

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strings"
)

// Seq 是一条带名称与描述的核酸序列（统一大写，不含空白）。
type Seq struct {
	Name string
	Desc string
	Data []byte
}

// ParseError 描述 FASTA 解析过程中发现的错误位置。
type ParseError struct {
	Line   int
	Column int
	Msg    string
}

func (e *ParseError) Error() string {
	loc := fmt.Sprintf("第 %d 行", e.Line)
	if e.Column > 0 {
		loc = fmt.Sprintf("第 %d 行第 %d 列", e.Line, e.Column)
	}
	return fmt.Sprintf("FASTA 解析错误（%s）：%s", loc, e.Msg)
}

// ParseFASTA 解析 FASTA 序列集合，非法字符或格式错误返回带位置的 *ParseError。
//
// 规则：
//   - 以 '>' 起始的行为序列头，首个空白前为名称、其余为描述；
//   - 以 ';' 起始的行（或首个非空白字符为 ';'）为注释，整行忽略；
//   - 空行忽略；序列行允许空白分隔，允许大小写混用，内部统一转为大写；
//   - 仅接受碱基 A/C/G/T（大小写不敏感），其余字符为非法字符；
//   - 序列数据前必须有头记录，且每条记录至少含一个碱基。
func ParseFASTA(r io.Reader) ([]Seq, error) {
	var seqs []Seq
	var cur *Seq
	var data bytes.Buffer
	headerLine := 0

	flush := func() error {
		if cur == nil {
			return nil
		}
		if data.Len() == 0 {
			return &ParseError{Line: headerLine, Msg: "序列记录缺少序列数据"}
		}
		cur.Data = append([]byte(nil), data.Bytes()...)
		seqs = append(seqs, *cur)
		cur = nil
		data.Reset()
		return nil
	}

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		raw := scanner.Text()
		trimmed := strings.TrimRight(raw, "\r\n")
		stripped := strings.TrimLeft(trimmed, " \t")
		prefix := len(trimmed) - len(stripped)

		if stripped == "" {
			continue
		}
		switch stripped[0] {
		case ';':
			continue
		case '>':
			if err := flush(); err != nil {
				return nil, err
			}
			header := strings.TrimSpace(stripped[1:])
			name := header
			desc := ""
			if i := strings.IndexAny(header, " \t"); i >= 0 {
				name = header[:i]
				desc = strings.TrimSpace(header[i+1:])
			}
			if name == "" {
				return nil, &ParseError{Line: lineNo, Column: 1, Msg: "序列头缺少名称"}
			}
			cur = &Seq{Name: name, Desc: desc}
			headerLine = lineNo
			continue
		}

		if cur == nil {
			return nil, &ParseError{Line: lineNo, Column: prefix + 1, Msg: "序列数据出现在任意序列头之前"}
		}
		for i := 0; i < len(trimmed); i++ {
			c := trimmed[i]
			if c == ' ' || c == '\t' {
				continue
			}
			up := c
			if up >= 'a' && up <= 'z' {
				up -= 'a' - 'A'
			}
			switch up {
			case 'A', 'C', 'G', 'T':
				data.WriteByte(up)
			default:
				return nil, &ParseError{
					Line:   lineNo,
					Column: i + 1,
					Msg:    fmt.Sprintf("非法字符 %q：仅允许碱基 A/C/G/T（不区分大小写）", c),
				}
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("读取 FASTA 失败：%w", err)
	}
	if err := flush(); err != nil {
		return nil, err
	}
	if len(seqs) == 0 {
		return nil, &ParseError{Msg: "输入为空或不含任何 FASTA 序列记录"}
	}
	return seqs, nil
}

// ParseFASTAStr 是 ParseFASTA 的字符串便捷封装。
func ParseFASTAStr(s string) ([]Seq, error) {
	return ParseFASTA(strings.NewReader(s))
}
