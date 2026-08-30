#!/usr/bin/env bash
set -euo pipefail

# Site administrators may override these two locations without modifying a
# release. The workflow root is resolved once so a running job remains pinned
# if the deployment's `current` symlink changes later.
WORKFLOW_ROOT="${CUTNRUN2TRACKS_ROOT:-/opt/bioinformatics/workflows/cutnrun2tracks/current}"
MAIN_ENV="${CUTNRUN2TRACKS_MAIN_ENV:-/opt/miniconda/envs/cutnrun2tracks-0.3.0}"
RESOLVED_WORKFLOW_ROOT="$(readlink -f -- "$WORKFLOW_ROOT" 2>/dev/null || true)"

[[ -n "$RESOLVED_WORKFLOW_ROOT" ]] || {
    echo "ERROR: cutnrun2tracks release cannot be resolved: $WORKFLOW_ROOT" >&2
    exit 127
}
WORKFLOW_ROOT="$RESOLVED_WORKFLOW_ROOT"
[[ -x "$WORKFLOW_ROOT/cutnrun2tracks.sh" ]] || {
    echo "ERROR: workflow is not executable: $WORKFLOW_ROOT/cutnrun2tracks.sh" >&2
    exit 127
}
[[ -x "$MAIN_ENV/bin/bash" && -x "$MAIN_ENV/bin/python3" ]] || {
    echo "ERROR: main environment is incomplete: $MAIN_ENV" >&2
    exit 127
}

# No `conda activate` is required and the caller's interactive shell is not
# modified. Optional versioned sidecars (SEACR, epic2, phantompeak/preseq) are
# exposed through administrator-managed launchers in /usr/local/bin.
export PATH="$MAIN_ENV/bin:/usr/local/bin:/usr/bin:/bin"
exec "$MAIN_ENV/bin/bash" "$WORKFLOW_ROOT/cutnrun2tracks.sh" "$@"
