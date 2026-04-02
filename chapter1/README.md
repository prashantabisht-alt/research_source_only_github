# chapter1

This folder is the cleanest Brownian / underdamped Langevin branch in the source-only repo.

## Best Entry Points

| File | Role |
| --- | --- |
| `SimulationFullLangevin2pointo_joint_enhanced50.f90` | best representative free-Langevin pipeline at `T_max = 50` |
| `SimulationFullLangevin2pointo_joint_enhanced5.f90` | same family at shorter time |
| `SimulationFullLangevin2pointo_joint_enhanced100.f90` | same family at longer time |
| `LangevinWithPotential.f90` | harmonic potential verification branch |
| `gaussian_simulation1.f90` | simple overdamped Gaussian benchmark |

## How To Run

```bash
bash run.sh --help
bash run.sh free50
bash run.sh compare50
```

## Notes

- Many similarly named files differ mainly by final time or overlay style.
- For one clean Brownian / Langevin reference branch, start with `SimulationFullLangevin2pointo_joint_enhanced50.f90`.
