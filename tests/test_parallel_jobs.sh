#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../scripts/lib/parallel_jobs.sh
source "$ROOT/scripts/lib/parallel_jobs.sh"

worker() {
    sleep "$1"
}

start_ns="$(date +%s%N)"
parallel_pool_init 2
parallel_pool_submit slow worker 0.8
parallel_pool_submit fast_first worker 0.1
parallel_pool_submit fast_second worker 0.6
parallel_pool_wait_all
end_ns="$(date +%s%N)"
elapsed_ms=$(( (end_ns - start_ns) / 1000000 ))

# A work-conserving pool completes in about 0.8 s. Waiting for the first
# submitted worker before refilling takes about 1.4 s.
(( elapsed_ms < 1200 )) || {
    echo "ERROR: parallel pool showed head-of-line blocking (${elapsed_ms} ms)" >&2
    exit 1
}

echo "Parallel pool scheduling test passed (${elapsed_ms} ms)"

failing_worker() {
    return 7
}

parallel_pool_init 2
parallel_pool_submit successful worker 0.1
parallel_pool_submit expected_failure failing_worker
if parallel_pool_wait_all; then
    echo "ERROR: parallel worker failure was not propagated" >&2
    exit 1
fi
[[ "$PARALLEL_POOL_FAILURES" -eq 1 ]]
[[ "${PARALLEL_POOL_FAILED_LABELS[*]}" == "expected_failure" ]]
echo "Parallel pool failure propagation test passed"
