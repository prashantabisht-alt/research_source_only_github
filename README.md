# TIFR DP1 Research Source-Only Repo

This repository is a source-only mirror of the larger DP1 workspace. It is meant to be the cleanest GitHub-facing version of the project for reading, collaboration, and future DP2 development.

Included here:

- Fortran source files
- gnuplot plotting scripts
- lightweight documentation added to make the repo easier to understand and reuse

Not included here:

- large generated `.txt` outputs
- PDFs and PNG figure exports
- report attachments and archive clutter from the full workspace

## Canonical Folders

These are the best folders to use first.

| Folder | Status | Main topic |
| --- | --- | --- |
| `chapter1` | canonical | underdamped Langevin / Brownian work |
| `ABP in 2D` | canonical | active Brownian particle results in 2D |
| `jerky particle` | canonical | jerky active particle / AOUP work |
| `RTP` | canonical | run-and-tumble particle work |
| `TCRW` | canonical | topological chiral random walker Python branch |

## Secondary Or Archive Folders

| Folder | Status | Notes |
| --- | --- | --- |
| `Brownian Motion` | early / exploratory | older Brownian and Langevin experiments |
| `Brownian motion separate` | exploratory | split variants and intermediate work |
| `Brownian Motion 2` | archive / duplicate-heavy | large archive with many repeated variants |
| `new_work` | later spillover | later CTRW / DSDT-oriented work |

## Added Cleanup Layer

This repo now includes:

- topic-level `README.md` files in the main folders
- portable `run.sh` wrappers for representative pipelines
- `FIGURE_MAP.md` for source -> plot mapping
- `ORGANIZATION_PLAN.md` for future cleanup guidance

## Quick Start

Recommended tools:

- `gfortran`
- `gnuplot`
- `python3`

The wrapper scripts compile with `-fno-range-check` because the bundled `mt.f90` RNG uses an older integer constant style.

Examples:

```bash
cd "chapter1"
bash run.sh free50

cd "../ABP in 2D"
bash run.sh fig1a

cd "../jerky particle"
bash run.sh fig3a

cd "../RTP"
bash run.sh fig1
```

For help:

```bash
bash run.sh --help
```

## Notes

- This is the best repo to share with Claude, Codex, or collaborators when you want code only.
- The full local workspace still contains the report, generated data, and additional archive material.
- If your main DP2 direction is jerky active particles, start with `jerky particle`, then use `chapter1` for Brownian / Langevin reference points.
- If you want the newer lattice-based chiral walker branch, start with `TCRW`.
