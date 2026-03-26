set logscale xy
set xlabel "t/τ_p"
set ylabel "MSD / (γ^2 v0^2 τ_p^6 / λ^2)"
set grid back lw 1
set border lw 1.5
set tics out
set key top left Left reverse samplen 2 spacing 1.1

# Styles
set style line 1 lc rgb "#1f77b4" pt 6 ps 1.3 lw 1.2    # sim (open circles)
set style line 2 lc rgb "#d62728" lw 2.4                # theory (solid line)
set style line 3 lc rgb "#2ca02c" lw 2 dt 2             # t^6 guide (dashed)
set style line 4 lc rgb "#9467bd" lw 2 dt 3             # t^5 guide (dash-dot)

plot \
  "msd_persistent_norm.txt"    u 1:2 w lp ls 1 t "simulation (OU)", \
  "theory_persistent_norm.txt" u 1:2 w  l  ls 2 t "theory (Eq. 23)", \
  (1.0/36.0)*x**6              w  l  ls 3 t "t^6 guide", \
  0.1*x**5                     w  l  ls 4 t "t^5 guide"
