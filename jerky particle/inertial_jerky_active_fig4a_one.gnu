set logscale xy
set xlabel "t / τI"
set ylabel "MSD / (γ^2 v0^2 τI^6 / λ^2)"
set grid back lw 1
set border lw 1.5
set tics out
set key top left spacing 1.2

# Styles
set style line 1 lc rgb "red"    pt 5 ps 0.9 lw 1.2    # simulation = red squares
set style line 2 lc rgb "black"  lw 3                  # theory = bold black line
set style line 3 lc rgb "#1f77b4" dt 2 lw 2            # t^6 guide = blue dashed
set style line 4 lc rgb "#2ca02c" dt 3 lw 2            # t^5 guide = green dash-dot
set style line 5 lc rgb "#9467bd" dt 4 lw 2            # t^3 guide = purple dotted

plot \
  "msd_active_taupratio_1e-3.txt"                u 1:2 w lp ls 1 t "simulation τp=τI/1000", \
  "msd_theory_taupratio_1e-3_overlay_numeric.txt" u 1:2 w  l  ls 2 t "theory τp=τI/1000", \
  0.0278*x**6    w l ls 3 t "t^6 prefactor", \
  0.0002*x**5    w l ls 4 t "t^5 prefactor", \
  0.001333*x**3  w l ls 5 t "t^3 prefactor"
