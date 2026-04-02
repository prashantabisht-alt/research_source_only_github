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
  fig1a         Compile and run the anisotropic MSD pipeline
  plot-fig1a    Open the Fig. 1(a) gnuplot script
  fig1b         Compile and run the snapshot pipeline
  plot-fig1b    Open the Fig. 1(b) gnuplot script
  fig2a         Compile and run the x/y sample pipeline
  plot-fig2a-x  Open the x-sample plot
  plot-fig2a-y  Open the y-sample plot
  fig2b         Compile and run the radial-distribution pipeline
  plot-fig2b    Open the radial-distribution plot
EOF
}

case "${1:-help}" in
  help|-h|--help)
    show_help
    ;;
  fig1a)
    compile_and_run "abp_fig1a_prash.f90" "abp_fig1a"
    ;;
  plot-fig1a)
    plot_gnu "abp_fig1a_prash.gnu"
    ;;
  fig1b)
    compile_and_run "abp_fig1b_snapshot.f90" "abp_fig1b_snapshot"
    ;;
  plot-fig1b)
    plot_gnu "abp_fig1b_snapshot.gnu"
    ;;
  fig2a)
    compile_and_run "abp_fig2a_samples.f90" "abp_fig2a_samples"
    ;;
  plot-fig2a-x)
    plot_gnu "abp_fig2a_samples_x.gnu"
    ;;
  plot-fig2a-y)
    plot_gnu "abp_fig2a_samples_y.gnu"
    ;;
  fig2b)
    compile_and_run "abp_fig2b_radial.f90" "abp_fig2b_radial"
    ;;
  plot-fig2b)
    plot_gnu "abp_fig2b_radial.gnu"
    ;;
  *)
    echo "Unknown target: $1" >&2
    show_help
    exit 1
    ;;
esac
