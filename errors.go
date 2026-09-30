package motif

import "fmt"

type taskError struct{ msg string }

func (e *taskError) Error() string { return e.msg }

func errUnknownTask(id string) error {
	return &taskError{msg: fmt.Sprintf("motif: 未知任务 %q", id)}
}

func errTaskExists(id string) error {
	return &taskError{msg: fmt.Sprintf("motif: 任务 %q 已存在（相同种子与配置）", id)}
}

func errBadStatus(id string, status TaskStatus) error {
	return &taskError{msg: fmt.Sprintf("motif: 任务 %q 当前状态 %s 不允许该操作", id, status)}
}
