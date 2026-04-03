#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

show_help() {
  cat <<'EOF'
Usage: bash run.sh <target>

Targets:
  core-check   Run the built-in sanity checks in tcrw_core.py
  phase1       Run the phase-1 PBC figure workflow
EOF
}

case "${1:-help}" in
  help|-h|--help)
    show_help
    ;;
  core-check)
    need_cmd python3
    (
      cd "$script_dir"
      python3 "tcrw_core.py"
    )
    ;;
  phase1)
    need_cmd python3
    (
      cd "$script_dir"
      python3 "tcrw_phase1_pbc.py"
    )
    ;;
  *)
    echo "Unknown target: $1" >&2
    show_help
    exit 1
    ;;
esac
