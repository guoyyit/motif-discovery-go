package motif

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// writeFixtureModule 在临时目录生成一个无外部依赖的迷你 Go 模块，
// fail 为 true 时其中的测试必然失败。
func writeFixtureModule(t *testing.T, fail bool) string {
	t.Helper()
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "go.mod"),
		[]byte("module fixture\n\ngo 1.21\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	body := `package fixture

import "testing"

func TestFixture(t *testing.T) {
`
	if fail {
		body += "\tt.Fatal(\"intentional failure\")\n"
	}
	body += "}\n"
	if err := os.WriteFile(filepath.Join(dir, "fixture_test.go"),
		[]byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	return dir
}

func runVerifyScript(t *testing.T, dir string) (string, int) {
	t.Helper()
	if _, err := exec.LookPath("sh"); err != nil {
		t.Skip("sh not available")
	}
	script, err := filepath.Abs("verify.sh")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(script); err != nil {
		t.Skip("verify.sh not found")
	}
	cmd := exec.Command("sh", script, dir)
	out, runErr := cmd.CombinedOutput()
	code := 0
	if runErr != nil {
		exitErr, ok := runErr.(*exec.ExitError)
		if !ok {
			t.Fatalf("cannot run verify.sh: %v", runErr)
		}
		code = exitErr.ExitCode()
	}
	return string(out), code
}

// TestVerifyScriptFailsOnFailingTests 验证流程在存在失败用例时：
// 退出码必须非零，且输出中不得出现任何表示全部通过的回显。
func TestVerifyScriptFailsOnFailingTests(t *testing.T) {
	out, code := runVerifyScript(t, writeFixtureModule(t, true))
	if code == 0 {
		t.Fatalf("verify.sh exited 0 despite failing tests; output:\n%s", out)
	}
	if strings.Contains(out, "VERIFY OK") {
		t.Fatalf("success banner printed despite failures; output:\n%s", out)
	}
	if !strings.Contains(out, "VERIFY FAILED") {
		t.Fatalf("failure summary missing; output:\n%s", out)
	}
}

// TestVerifyScriptSucceedsOnPassingTests 全部通过时才允许输出成功结论，
// 且退出码必须为 0（汇总信息与退出码一致）。
func TestVerifyScriptSucceedsOnPassingTests(t *testing.T) {
	out, code := runVerifyScript(t, writeFixtureModule(t, false))
	if code != 0 {
		t.Fatalf("verify.sh exited %d on passing module; output:\n%s", code, out)
	}
	if !strings.Contains(out, "VERIFY OK") {
		t.Fatalf("success banner missing on passing module; output:\n%s", out)
	}
}
