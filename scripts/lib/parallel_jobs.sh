#!/usr/bin/env bash
set -euo pipefail

parallel_require_positive_integer() {
    local name="${1:?name required}" value="${2:-}"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: $name must be a positive integer; found ${value:-<empty>}" >&2
        return 1
    }
}

parallel_pool_init() {
    parallel_require_positive_integer "parallel job limit" "${1:?limit required}"
    PARALLEL_POOL_MAX="$1"
    PARALLEL_POOL_FAILURES=0
    PARALLEL_POOL_PIDS=()
    PARALLEL_POOL_LABELS=()
    PARALLEL_POOL_FAILED_LABELS=()
}

parallel_pool_wait_one() {
    (( ${#PARALLEL_POOL_PIDS[@]} > 0 )) || return 0
    local finished_pid="" status=0 index label=""
    if wait -n -p finished_pid "${PARALLEL_POOL_PIDS[@]}"; then
        status=0
    else
        status=$?
    fi
    [[ -n "$finished_pid" ]] || {
        echo "ERROR: parallel pool could not identify the completed worker" >&2
        return 1
    }
    for index in "${!PARALLEL_POOL_PIDS[@]}"; do
        if [[ "${PARALLEL_POOL_PIDS[$index]}" == "$finished_pid" ]]; then
            label="${PARALLEL_POOL_LABELS[$index]}"
            unset "PARALLEL_POOL_PIDS[$index]" "PARALLEL_POOL_LABELS[$index]"
            PARALLEL_POOL_PIDS=("${PARALLEL_POOL_PIDS[@]}")
            PARALLEL_POOL_LABELS=("${PARALLEL_POOL_LABELS[@]}")
            break
        fi
    done
    [[ -n "$label" ]] || {
        echo "ERROR: completed worker $finished_pid is absent from the parallel pool" >&2
        return 1
    }
    if (( status != 0 )); then
        PARALLEL_POOL_FAILURES=$((PARALLEL_POOL_FAILURES + 1))
        PARALLEL_POOL_FAILED_LABELS+=("$label")
    fi
}

parallel_pool_submit() {
    local label="${1:?label required}"
    shift
    (( $# > 0 )) || return 1
    while (( ${#PARALLEL_POOL_PIDS[@]} >= PARALLEL_POOL_MAX )); do
        parallel_pool_wait_one
    done
    "$@" &
    PARALLEL_POOL_PIDS+=("$!")
    PARALLEL_POOL_LABELS+=("$label")
}

parallel_pool_wait_all() {
    while (( ${#PARALLEL_POOL_PIDS[@]} > 0 )); do
        parallel_pool_wait_one
    done
    if (( PARALLEL_POOL_FAILURES > 0 )); then
        local joined
        joined="$(IFS=,; echo "${PARALLEL_POOL_FAILED_LABELS[*]}")"
        echo "ERROR: parallel workers failed: $joined" >&2
        return 1
    fi
}
