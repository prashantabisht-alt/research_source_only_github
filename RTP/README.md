# RTP

This folder contains the cleanest run-and-tumble particle branch in the source-only repo.

## Main Files

| File | Role |
| --- | --- |
| `rtp_fig1.f90` | Monte Carlo generation of `P(x,t)` for multiple `D` values |
| `rtp_fig1_theory.py` | theory data generation for the same panels |
| `rtp_fig1_Dpoint03.gnu` | plotting script for `D = 0.03` |
| `rtp_fig1_Dpoint1.gnu` | plotting script for `D = 0.10` |
| `rtp_fig1_Dpoint2.gnu` | plotting script for `D = 0.20` |

## How To Run

```bash
bash run.sh --help
bash run.sh fig1
bash run.sh plot-d003
```

## Notes

- This is a compact simulation-plus-theory branch and a good example of a clean figure pipeline.
