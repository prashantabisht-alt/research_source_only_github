# Figure And Pipeline Map

This file records representative source -> plot chains for the main code branches in the source-only repo.

## chapter1

| Topic | Source | Plot layer | Notes |
| --- | --- | --- | --- |
| Free underdamped Langevin at `T_max = 50` | `chapter1/SimulationFullLangevin2pointo_joint_enhanced50.f90` | `chapter1/SimulationFullLangevin2pointo_joint_enhanced50.gnu`, `chapter1/SimulationFullLangevin2pointo_joint_enhanced50_compare.gnu`, `chapter1/SimulationFullLangevin2pointo_joint_enhanced_compare_theory.gnu` | Best representative Brownian/Langevin branch |
| Same family at other times | `chapter1/SimulationFullLangevin2pointo_joint_enhanced5.f90`, `chapter1/SimulationFullLangevin2pointo_joint_enhanced100.f90` | matching `*.gnu` files | Same pipeline at different final times |
| Langevin with harmonic potential | `chapter1/LangevinWithPotential.f90` | nearby potential-related scripts | Good equilibrium-check branch |
| Overdamped Gaussian benchmark | `chapter1/gaussian_simulation1.f90` | `chapter1/overdamped_pdf.gnu` | Simple diffusion benchmark |

## ABP in 2D

| Topic | Source | Plot layer | Notes |
| --- | --- | --- | --- |
| Fig. 1(a): anisotropic MSD | `ABP in 2D/abp_fig1a_prash.f90` | `ABP in 2D/abp_fig1a_prash.gnu` | Clean figure-driven branch |
| Fig. 1(b): snapshot | `ABP in 2D/abp_fig1b_snapshot.f90` | `ABP in 2D/abp_fig1b_snapshot.gnu`, `ABP in 2D/abp_1b_snapshot_raw.gnu`, `ABP in 2D/abp_1b_snapshot_3d.gnu` | Spatial visualization pipeline |
| Fig. 2(a): x/y samples | `ABP in 2D/abp_fig2a_samples.f90` | `ABP in 2D/abp_fig2a_samples_x.gnu`, `ABP in 2D/abp_fig2a_samples_y.gnu` | Distribution slices |
| Fig. 2(b): radial distribution | `ABP in 2D/abp_fig2b_radial.f90` | `ABP in 2D/abp_fig2b_radial.gnu` | Radial theory/simulation branch |

## jerky particle

| Topic | Source | Plot layer | Notes |
| --- | --- | --- | --- |
| Fig. 3(a): pure jerky AOUP MSD | `jerky particle/aoup_fig3a.f90` | `jerky particle/aoup_fig3a.gnu` | Best single-particle jerky entry point |
| Fig. 3(b): exponent from MSD | `jerky particle/aoup_fig3b.f90` | `jerky particle/aoup_fig3b.gnu` | Post-processing branch |
| Fig. 4(a): fast-noise branch | `jerky particle/inertial_jerky_active_fig4a_one.f90` | `jerky particle/inertial_jerky_active_fig4a_one.gnu` | Fast-noise crossover |
| Fig. 4(a): slow-noise branch | `jerky particle/inertial_jerky_active_fig4a_two.f90` | `jerky particle/inertial_jerky_active_fig4a_two.gnu` | Slow-noise crossover |
| Fig. 4(b): exponent analysis | `jerky particle/aoup_fig4b.f90` | `jerky particle/aoup_fig4b_one.gnu`, `jerky particle/aoup_fig4b_two.gnu` | Derived exponent plots |
| Eq. 6-8 diagnostics | `jerky particle/aoup_eq678.f90` | `jerky particle/aoup_eq678_6.gnu`, `jerky particle/aoup_eq678_7.gnu` | Kernel / correlation diagnostics |

## RTP

| Topic | Source | Plot layer | Notes |
| --- | --- | --- | --- |
| Fig. 1 Monte Carlo | `RTP/rtp_fig1.f90` | `RTP/rtp_fig1_Dpoint03.gnu`, `RTP/rtp_fig1_Dpoint1.gnu`, `RTP/rtp_fig1_Dpoint2.gnu` | Simulation branch |
| Fig. 1 theory overlay | `RTP/rtp_fig1_theory.py` | same gnuplot files | Theory and MC share the same plotting scripts |
