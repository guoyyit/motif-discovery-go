#!/bin/sh
# 本地验收流程：构建、静态检查与测试全部通过才输出成功结论。
# 任一步骤失败即以非零退出码结束，且不会打印任何表示全部通过的回显。
#
# 用法：sh verify.sh [目录]
#   目录缺省为脚本所在仓库根目录；测试可传入其他模块目录以验证脚本行为。
set -u

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TARGET=${1:-"$ROOT"}

if ! cd "$TARGET"; then
    echo "verify: cannot enter directory: $TARGET" >&2
    exit 2
fi

failed=0

run_step() {
    step_name=$1
    shift
    echo "== ${step_name}"
    if "$@"; then
        echo "-- ${step_name}: OK"
    else
        rc=$?
        echo "-- ${step_name}: FAILED (exit ${rc})" >&2
        failed=1
    fi
}

run_step "build" go build ./...
run_step "vet" go vet ./...
run_step "test" go test -count=1 ./...

if [ "${failed}" -ne 0 ]; then
    echo "VERIFY FAILED: one or more steps did not pass" >&2
    exit 1
fi

echo "VERIFY OK: build, vet and tests all passed"
exit 0
