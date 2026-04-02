# jerky particle

This folder is the strongest jerky-active-particle branch in the source-only repo and the best starting point for jerky-focused DP2 work.

## Best Entry Points

| File | Role |
| --- | --- |
| `aoup_fig3a.f90` | clean Fig. 3(a) MSD pipeline with theory output |
| `aoup_fig3b.f90` | post-processing of Fig. 3(a) into an effective exponent |
| `inertial_jerky_active_fig4a_one.f90` | fast-noise branch for Fig. 4(a) |
| `inertial_jerky_active_fig4a_two.f90` | slow-noise branch for Fig. 4(a) |
| `aoup_fig4b.f90` | exponent post-processing for the Fig. 4(a) branches |
| `aoup_eq678.f90` | kernel / correlation diagnostics |

## How To Run

```bash
bash run.sh --help
bash run.sh fig3a
bash run.sh all-fig4
```

## Notes

- If your DP2 is jerky-focused, this folder should be your primary base.
