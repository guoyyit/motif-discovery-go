#!/usr/bin/env bash
# 本地验收入口（无外部服务、无网络）。
#
# 语义约定（必须严格保持）：
#   - 构建、vet、任一测试用例失败时：立即以非零退出码结束，
#     除 go 工具自身输出的失败信息外，绝不打印任何“全部通过/成功”回显；
#   - 仅当全部检查真正通过时，才打印唯一的成功标记 ALL_CHECKS_PASSED
#     并以退出码 0 结束。
#
# 用法：scripts/verify.sh [模块目录]（默认为脚本所在仓库根目录）。
set -u
set -o pipefail

DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$DIR" || {
	echo "verify: cannot enter module directory: $DIR" >&2
	exit 2
}

# GOPROXY=off 与 CGO 默认值保证离线运行；不访问任何外部服务。
export GOPROXY=off
export GOFLAGS=-mod=mod
export GOTOOLCHAIN=local

run() {
	echo "==> $*"
	if ! "$@"; then
		echo "verify: FAILED: $*" >&2
		exit 1
	fi
}

run go build ./...
run go vet ./...
run go test -count=1 ./...

# 只有上面全部成功才会执行到这里。
echo "ALL_CHECKS_PASSED"
