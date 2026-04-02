#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="$script_dir/.build"
mkdir -p "$build_dir"

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

compile_and_run() {
  local src="$1"
  local exe="$2"
  need_cmd gfortran
  (
    cd "$script_dir"
    gfortran -O2 -fno-range-check "$src" -o "$build_dir/$exe"
    "$build_dir/$exe"
  )
}

plot_gnu() {
  local script="$1"
  need_cmd gnuplot
  (
    cd "$script_dir"
    gnuplot "$script"
  )
}

show_help() {
  cat <<'EOF'
Usage: bash run.sh <target>

Targets:
  fig1         Compile and run the RTP Monte Carlo generator, then build theory files
  fig1-sim     Compile and run the RTP Monte Carlo generator only
  fig1-theory  Run the theory generator only
  plot-d003    Open the D=0.03 plot
  plot-d010    Open the D=0.10 plot
  plot-d020    Open the D=0.20 plot
  all-plots    Open all three RTP figure plots in sequence
EOF
}

case "${1:-help}" in
  help|-h|--help)
    show_help
    ;;
  fig1)
    compile_and_run "rtp_fig1.f90" "rtp_fig1"
    need_cmd python3
    (
      cd "$script_dir"
      python3 "rtp_fig1_theory.py"
    )
    ;;
  fig1-sim)
    compile_and_run "rtp_fig1.f90" "rtp_fig1"
    ;;
  fig1-theory)
    need_cmd python3
    (
      cd "$script_dir"
      python3 "rtp_fig1_theory.py"
    )
    ;;
  plot-d003)
    plot_gnu "rtp_fig1_Dpoint03.gnu"
    ;;
  plot-d010)
    plot_gnu "rtp_fig1_Dpoint1.gnu"
    ;;
  plot-d020)
    plot_gnu "rtp_fig1_Dpoint2.gnu"
    ;;
  all-plots)
    plot_gnu "rtp_fig1_Dpoint03.gnu"
    plot_gnu "rtp_fig1_Dpoint1.gnu"
    plot_gnu "rtp_fig1_Dpoint2.gnu"
    ;;
  *)
    echo "Unknown target: $1" >&2
    show_help
    exit 1
    ;;
esac
