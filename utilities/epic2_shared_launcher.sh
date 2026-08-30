#!/usr/bin/env bash
set -euo pipefail

EPIC2_ENV="${CUTNRUN2TRACKS_EPIC2_ENV:-/opt/miniconda/envs/cutnrun2tracks-epic2-0.3.0}"
[[ -x "$EPIC2_ENV/bin/epic2" ]] || {
    echo "ERROR: epic2 sidecar is incomplete: $EPIC2_ENV" >&2
    exit 127
}
export PATH="$EPIC2_ENV/bin:/usr/bin:/bin"
exec "$EPIC2_ENV/bin/epic2" "$@"
