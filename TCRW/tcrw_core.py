"""
TCRW Core Simulation Engine
============================
Topological Chiral Random Walker on a 2D square lattice.

Model rules (Osat, Meyberg, Metson & Speck, arXiv:2602.12020):
  State: (x, y, d) where d in {0,1,2,3} = {up, right, down, left}  (CCW ordering)

  At each discrete time step:
    With prob D_r:  NOISE STEP
      - Walker stays at (x,y)
      - Director rotates: CCW (d -> d+1 mod 4) with prob omega
                          CW  (d -> d-1 mod 4) with prob (1-omega)

    With prob (1-D_r):  CHIRAL STEP
      - Walker translates one step in direction d
      - Then director rotates: CW  (d -> d-1 mod 4) with prob omega
                               CCW (d -> d+1 mod 4) with prob (1-omega)
      NOTE: rotation chirality in chiral step is opposite to noise step.

      OBC special rule: if translation would move walker off-grid,
      the entire chiral step is blocked (no translation and no rotation).

Author: Prashant Bisht, TIFR Hyderabad
"""

import numpy as np

DX = np.array([0, 1, 0, -1], dtype=np.int64)
DY = np.array([1, 0, -1, 0], dtype=np.int64)


def simulate_tcrw_pbc(omega, D_r, L, T_steps, N_traj=1, seed=42,
                      record_traj=False, record_interval=1):
    """
    Simulate TCRW on an LxL lattice with periodic boundary conditions.
    """
    rng = np.random.default_rng(seed)

    x = np.full(N_traj, L // 2, dtype=np.int64)
    y = np.full(N_traj, L // 2, dtype=np.int64)
    d = rng.integers(0, 4, size=N_traj)

    ux = np.zeros(N_traj, dtype=np.int64)
    uy = np.zeros(N_traj, dtype=np.int64)

    n_records = T_steps // record_interval + 1
    msd = np.zeros(n_records, dtype=np.float64)
    times = np.zeros(n_records, dtype=np.int64)

    if record_traj:
        traj_x = np.zeros(n_records, dtype=np.int64)
        traj_y = np.zeros(n_records, dtype=np.int64)
        traj_x[0] = ux[0]
        traj_y[0] = uy[0]

    rec_idx = 1

    for t in range(1, T_steps + 1):
        r_step = rng.random(N_traj)
        r_rot = rng.random(N_traj)

        is_noise = r_step < D_r
        is_chiral = ~is_noise

        noise_ccw = is_noise & (r_rot < omega)
        noise_cw = is_noise & (r_rot >= omega)
        d[noise_ccw] = (d[noise_ccw] + 1) % 4
        d[noise_cw] = (d[noise_cw] - 1) % 4

        step_dx = DX[d]
        step_dy = DY[d]

        x[is_chiral] = (x[is_chiral] + step_dx[is_chiral]) % L
        y[is_chiral] = (y[is_chiral] + step_dy[is_chiral]) % L
        ux[is_chiral] += step_dx[is_chiral]
        uy[is_chiral] += step_dy[is_chiral]

        chiral_cw = is_chiral & (r_rot < omega)
        chiral_ccw = is_chiral & (r_rot >= omega)
        d[chiral_cw] = (d[chiral_cw] - 1) % 4
        d[chiral_ccw] = (d[chiral_ccw] + 1) % 4

        if t % record_interval == 0:
            r2 = ux.astype(np.float64) ** 2 + uy.astype(np.float64) ** 2
            msd[rec_idx] = np.mean(r2)
            times[rec_idx] = t
            if record_traj:
                traj_x[rec_idx] = ux[0]
                traj_y[rec_idx] = uy[0]
            rec_idx += 1

    result = {
        "msd": msd[:rec_idx],
        "times": times[:rec_idx],
        "final_x": x,
        "final_y": y,
        "final_d": d,
    }
    if record_traj:
        result["traj_x"] = traj_x[:rec_idx]
        result["traj_y"] = traj_y[:rec_idx]
    return result


def simulate_tcrw_obc(omega, D_r, L, T_steps, N_traj=1, seed=42,
                      record_traj=False, record_interval=1,
                      track_visits=False, track_currents=False,
                      track_step_type=False):
    """
    Simulate TCRW on an LxL lattice with open (hard-wall) boundary conditions.
    """
    rng = np.random.default_rng(seed)

    x = rng.integers(0, L, size=N_traj)
    y = rng.integers(0, L, size=N_traj)
    d = rng.integers(0, 4, size=N_traj)

    x0 = x.copy()
    y0 = y.copy()

    if track_visits:
        visits = np.zeros((L, L), dtype=np.int64)

    if track_currents:
        Jx = np.zeros((L, L), dtype=np.float64)
        Jy = np.zeros((L, L), dtype=np.float64)
        Jx_Dr = np.zeros((L, L), dtype=np.float64)
        Jy_Dr = np.zeros((L, L), dtype=np.float64)
        Jx_omega = np.zeros((L, L), dtype=np.float64)
        Jy_omega = np.zeros((L, L), dtype=np.float64)
        prev_was_noise = np.zeros(N_traj, dtype=bool)

    n_records = T_steps // record_interval + 1
    msd = np.zeros(n_records, dtype=np.float64)
    times = np.zeros(n_records, dtype=np.int64)

    if record_traj:
        traj_x = np.zeros(n_records, dtype=np.int64)
        traj_y = np.zeros(n_records, dtype=np.int64)
        traj_x[0] = x[0]
        traj_y[0] = y[0]

    rec_idx = 1

    for t in range(1, T_steps + 1):
        r_step = rng.random(N_traj)
        r_rot = rng.random(N_traj)

        is_noise = r_step < D_r
        is_chiral = ~is_noise

        if track_visits:
            for i in range(N_traj):
                visits[x[i], y[i]] += 1

        noise_ccw = is_noise & (r_rot < omega)
        noise_cw = is_noise & (r_rot >= omega)
        d[noise_ccw] = (d[noise_ccw] + 1) % 4
        d[noise_cw] = (d[noise_cw] - 1) % 4

        if np.any(is_chiral):
            new_x = x + DX[d]
            new_y = y + DY[d]
            can_move = is_chiral & (new_x >= 0) & (new_x < L) & (new_y >= 0) & (new_y < L)

            if track_currents:
                for i in np.where(can_move)[0]:
                    ddx = DX[d[i]]
                    ddy = DY[d[i]]
                    Jx[x[i], y[i]] += ddx
                    Jy[x[i], y[i]] += ddy
                    if prev_was_noise[i]:
                        Jx_Dr[x[i], y[i]] += ddx
                        Jy_Dr[x[i], y[i]] += ddy
                    else:
                        Jx_omega[x[i], y[i]] += ddx
                        Jy_omega[x[i], y[i]] += ddy

            x[can_move] = new_x[can_move]
            y[can_move] = new_y[can_move]

            chiral_cw = can_move & (r_rot < omega)
            chiral_ccw = can_move & (r_rot >= omega)
            d[chiral_cw] = (d[chiral_cw] - 1) % 4
            d[chiral_ccw] = (d[chiral_ccw] + 1) % 4

        if track_currents:
            prev_was_noise[:] = False
            prev_was_noise[is_noise] = True

        if t % record_interval == 0:
            r2 = (x.astype(np.float64) - x0.astype(np.float64)) ** 2 + \
                 (y.astype(np.float64) - y0.astype(np.float64)) ** 2
            msd[rec_idx] = np.mean(r2)
            times[rec_idx] = t
            if record_traj:
                traj_x[rec_idx] = x[0]
                traj_y[rec_idx] = y[0]
            rec_idx += 1

    if track_visits:
        for i in range(N_traj):
            visits[x[i], y[i]] += 1

    result = {
        "msd": msd[:rec_idx],
        "times": times[:rec_idx],
        "final_x": x,
        "final_y": y,
        "final_d": d,
    }
    if record_traj:
        result["traj_x"] = traj_x[:rec_idx]
        result["traj_y"] = traj_y[:rec_idx]
    if track_visits:
        result["visits"] = visits
    if track_currents:
        total_steps = T_steps * N_traj
        result["Jx"] = Jx / total_steps
        result["Jy"] = Jy / total_steps
        result["Jx_Dr"] = Jx_Dr / total_steps
        result["Jy_Dr"] = Jy_Dr / total_steps
        result["Jx_omega"] = Jx_omega / total_steps
        result["Jy_omega"] = Jy_omega / total_steps
    return result


def measure_diffusion_coeff(omega, D_r, L=200, T_steps=500000, N_traj=500,
                            seed=42, fit_frac=0.5):
    """
    Measure the diffusion coefficient D from the long-time MSD slope.
    MSD = 4Dt in 2D.
    """
    res = simulate_tcrw_pbc(omega, D_r, L, T_steps, N_traj, seed,
                            record_interval=100)
    t = res["times"].astype(np.float64)
    msd = res["msd"]

    n = len(t)
    i0 = int(n * (1 - fit_frac))
    t_fit = t[i0:]
    msd_fit = msd[i0:]
    slope = np.polyfit(t_fit, msd_fit, 1)[0]
    return slope / 4.0


if __name__ == "__main__":
    print("=== TCRW Core Sanity Checks ===\n")

    D = measure_diffusion_coeff(omega=0.5, D_r=0.001, L=200, T_steps=200000,
                                N_traj=200, seed=42)
    print(f"omega=0.5, D_r=0.001: D = {D:.4f}")

    res = simulate_tcrw_pbc(omega=0.0, D_r=0.0, L=10, T_steps=8, N_traj=1,
                            seed=42, record_traj=True, record_interval=1)
    print("\nD_r=0, omega=0 (det. CCW chiral rotor):")
    print(f"  traj_x: {res['traj_x']}")
    print(f"  traj_y: {res['traj_y']}")
    print("  Should trace CCW orbit with period 4")

    res = simulate_tcrw_pbc(omega=1.0, D_r=1.0, L=10, T_steps=1000, N_traj=100,
                            seed=42)
    print(f"\nD_r=1, omega=1 (pure spinor): final MSD = {res['msd'][-1]:.6f}")
    print("  Should be 0 (no translation)")

    print("\nDone.")
