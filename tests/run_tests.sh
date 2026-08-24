#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "$ROOT/tests/check_bash_syntax.sh"
bash "$ROOT/tests/test_parallel_jobs.sh"
bash "$ROOT/tests/test_filter_cleanup.sh"
bash "$ROOT/tests/test_peakcall_tolerance.sh"
python3 -m unittest discover -s "$ROOT/tests" -p 'test_*.py' -v
echo "All cutnrun2tracks tests passed"
