#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapfile -t scripts < <(find "$ROOT" -type f -name '*.sh' -print | sort)
(( ${#scripts[@]} > 0 ))
for script in "${scripts[@]}"; do bash -n "$script"; done
echo "Bash syntax OK: ${#scripts[@]} files"
