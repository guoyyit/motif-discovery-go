package motif_test

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// scriptPath 返回仓库内 scripts/verify.sh 的绝对路径。
func scriptPath(t *testing.T) string {
	t.Helper()
	wd, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	return filepath.Join(wd, "scripts", "verify.sh")
}

// runVerify 在指定模块目录执行验收脚本，返回退出码与合并输出。
func runVerify(t *testing.T, moduleDir string) (int, string) {
	t.Helper()
	cmd := exec.Command("/bin/bash", scriptPath(t), moduleDir)
	out, _ := cmd.CombinedOutput()
	code := cmd.ProcessState.ExitCode()
	return code, string(out)
}

// TestVerifyScriptFailsNonZeroNoSuccessEcho 需求 7：存在失败用例时
// 脚本必须非零退出，且绝不在末尾打印任何表示全部通过的回显。
func TestVerifyScriptFailsNonZeroNoSuccessEcho(t *testing.T) {
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "go.mod"),
		[]byte("module failingexample\n\ngo 1.25.0\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(dir, "x_test.go"), []byte(
		"package failingexample\n\nimport \"testing\"\n\n"+
			"func TestAlwaysFails(t *testing.T) {\n\tt.Fatal(\"synthetic failure\")\n}\n"),
		0o644); err != nil {
		t.Fatal(err)
	}

	code, out := runVerify(t, dir)
	if code == 0 {
		t.Fatalf("verify script must exit non-zero when a test fails; output:\n%s", out)
	}
	if strings.Contains(out, "ALL_CHECKS_PASSED") {
		t.Fatalf("verify script must not print success echo on failure; output:\n%s", out)
	}
	if !strings.Contains(out, "FAIL") {
		t.Fatalf("failure output should report FAIL; output:\n%s", out)
	}
}

// TestVerifyScriptSucceedsOnPassingModule 全部通过时退出码为 0，
// 且输出唯一成功标记。
func TestVerifyScriptSucceedsOnPassingModule(t *testing.T) {
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "go.mod"),
		[]byte("module passingexample\n\ngo 1.25.0\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(dir, "x_test.go"), []byte(
		"package passingexample\n\nimport \"testing\"\n\n"+
			"func TestAlwaysPasses(t *testing.T) {\n\tif 1+1 != 2 {\n\t\tt.Fatal(\"math broken\")\n\t}\n}\n"),
		0o644); err != nil {
		t.Fatal(err)
	}

	code, out := runVerify(t, dir)
	if code != 0 {
		t.Fatalf("verify script should exit 0 on passing module (code=%d); output:\n%s",
			code, out)
	}
	if !strings.Contains(out, "ALL_CHECKS_PASSED") {
		t.Fatalf("passing run must contain success marker; output:\n%s", out)
	}
}

// TestVerifyScriptFailsOnBuildError 构建失败同样必须非零退出且无成功回显。
func TestVerifyScriptFailsOnBuildError(t *testing.T) {
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "go.mod"),
		[]byte("module buildbroken\n\ngo 1.25.0\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(dir, "x.go"),
		[]byte("package buildbroken\n\nfunc broken( {}\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	code, out := runVerify(t, dir)
	if code == 0 {
		t.Fatalf("build failure must yield non-zero exit; output:\n%s", out)
	}
	if strings.Contains(out, "ALL_CHECKS_PASSED") {
		t.Fatalf("build failure must not print success echo; output:\n%s", out)
	}
}
