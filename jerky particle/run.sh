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
  fig3a             Compile and run the clean Fig. 3(a) pipeline
  plot-fig3a        Open the Fig. 3(a) gnuplot script
  fig3b             Compile and run the alpha-from-MSD post-processing
  plot-fig3b        Open the Fig. 3(b) gnuplot script
  fig4a-fast        Compile and run the fast-noise Fig. 4(a) branch
  fig4a-slow        Compile and run the slow-noise Fig. 4(a) branch
  fig4b             Compile and run the Fig. 4(b) post-processing
  all-fig4          Run fast-noise, slow-noise, then alpha post-processing
  plot-fig4-fast    Open the fast-noise Fig. 4(a) plot
  plot-fig4-slow    Open the slow-noise Fig. 4(a) plot
  plot-fig4-alpha1  Open the fast-noise alpha plot
  plot-fig4-alpha2  Open the slow-noise alpha plot
  eq678             Compile and run the kernel / autocorrelation diagnostics
EOF
}

case "${1:-help}" in
  help|-h|--help)
    show_help
    ;;
  fig3a)
    compile_and_run "aoup_fig3a.f90" "jerky_fig3a"
    ;;
  plot-fig3a)
    plot_gnu "aoup_fig3a.gnu"
    ;;
  fig3b)
    compile_and_run "aoup_fig3b.f90" "jerky_fig3b"
    ;;
  plot-fig3b)
    plot_gnu "aoup_fig3b.gnu"
    ;;
  fig4a-fast)
    compile_and_run "inertial_jerky_active_fig4a_one.f90" "jerky_fig4a_fast"
    ;;
  fig4a-slow)
    compile_and_run "inertial_jerky_active_fig4a_two.f90" "jerky_fig4a_slow"
    ;;
  fig4b)
    compile_and_run "aoup_fig4b.f90" "jerky_fig4b"
    ;;
  all-fig4)
    compile_and_run "inertial_jerky_active_fig4a_one.f90" "jerky_fig4a_fast"
    compile_and_run "inertial_jerky_active_fig4a_two.f90" "jerky_fig4a_slow"
    compile_and_run "aoup_fig4b.f90" "jerky_fig4b"
    ;;
  plot-fig4-fast)
    plot_gnu "inertial_jerky_active_fig4a_one.gnu"
    ;;
  plot-fig4-slow)
    plot_gnu "inertial_jerky_active_fig4a_two.gnu"
    ;;
  plot-fig4-alpha1)
    plot_gnu "aoup_fig4b_one.gnu"
    ;;
  plot-fig4-alpha2)
    plot_gnu "aoup_fig4b_two.gnu"
    ;;
  eq678)
    compile_and_run "aoup_eq678.f90" "jerky_eq678"
    ;;
  *)
    echo "Unknown target: $1" >&2
    show_help
    exit 1
    ;;
esac
