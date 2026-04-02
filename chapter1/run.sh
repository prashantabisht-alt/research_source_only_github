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
  free50          Compile and run the main free-Langevin T=50 pipeline
  plot-free50     Open the basic joint-distribution plot
  compare50       Open the simulation-vs-theory heatmap comparison
  theory-panels   Open the multi-time theory heatmap panel
  potential       Compile and run the harmonic-potential verification code
  overdamped      Compile and run the simple Gaussian benchmark
  plot-overdamped Open the overdamped PDF plot
EOF
}

case "${1:-help}" in
  help|-h|--help)
    show_help
    ;;
  free50)
    compile_and_run "SimulationFullLangevin2pointo_joint_enhanced50.f90" "chapter1_free50"
    ;;
  plot-free50)
    plot_gnu "SimulationFullLangevin2pointo_joint_enhanced50.gnu"
    ;;
  compare50)
    plot_gnu "SimulationFullLangevin2pointo_joint_enhanced50_compare.gnu"
    ;;
  theory-panels)
    plot_gnu "SimulationFullLangevin2pointo_joint_enhanced_compare_theory.gnu"
    ;;
  potential)
    compile_and_run "LangevinWithPotential.f90" "chapter1_potential"
    ;;
  overdamped)
    compile_and_run "gaussian_simulation1.f90" "chapter1_overdamped"
    ;;
  plot-overdamped)
    plot_gnu "overdamped_pdf.gnu"
    ;;
  *)
    echo "Unknown target: $1" >&2
    show_help
    exit 1
    ;;
esac
