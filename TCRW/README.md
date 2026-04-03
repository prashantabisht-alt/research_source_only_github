# TCRW

This folder contains a Python implementation of a topological chiral random walker (TCRW) on a square lattice.

## Files

| File | Role |
| --- | --- |
| `tcrw_core.py` | main simulation engine with PBC and OBC variants |
| `tcrw_phase1_pbc.py` | figure-style driver for periodic-boundary-condition results |
| `run.sh` | simple wrapper to run the Python checks and phase-1 script |

## Requirements

- `python3`
- `numpy`
- `matplotlib`

## How To Run

```bash
bash run.sh --help
bash run.sh core-check
bash run.sh phase1
```

## Notes

- The local full workspace also contains generated PNG figures from this branch, but they are intentionally not included in the source-only repo.
- This branch is separate from the older Fortran/gnuplot pipelines and is one of the few Python-first parts of the project.
