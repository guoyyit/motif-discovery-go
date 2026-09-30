package motif

import (
	"context"
	"testing"
	"time"
)

func TestTaskCancel(t *testing.T) {
	rng := NewRNG(7171)
	seqs := generateSequences(rng, 20, 100)

	cfg := testConfig(1000, 1000) // 计算量足够大，确保取消时仍在运行
	mgr := NewTaskManager()
	id, err := mgr.Submit(seqs, cfg)
	if err != nil {
		t.Fatal(err)
	}
	time.Sleep(5 * time.Millisecond)
	if err := mgr.Cancel(context.Background(), id); err != nil {
		t.Fatal(err)
	}

	deadline := time.Now().Add(10 * time.Second)
	var status TaskStatus
	for time.Now().Before(deadline) {
		status, _, _ = mgr.Get(id)
		if status == StatusCanceled || status == StatusFailed {
			break
		}
		time.Sleep(2 * time.Millisecond)
	}
	if status != StatusCanceled {
		t.Fatalf("expected canceled, got %s", status)
	}
}
