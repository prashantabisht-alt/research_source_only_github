# ===== Velocity & Position PDFs side-by-side =====
reset
set term qt size 1600,650
set multiplot layout 1,2 margins 0.07,0.97,0.12,0.93 spacing 0.06,0.02 \
    title "Langevin (inertial, V=0): Velocity & Position PDFs (T_max = 50)"

# Common look
set grid back lw 1
set border lw 1.6
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 600
set style line 1 lc rgb "#d62728" pt 5 ps 0.9 lw 1.4   # vel sim
set style line 2 lc rgb "black"   lw 2.8               # vel theory
set style line 3 lc rgb "#1f77b4" pt 7 ps 0.9 lw 1.4   # pos sim
set style line 4 lc rgb "black"   lw 2.8               # pos theory

# Left panel: Velocity PDF
set title "Velocity PDF (Maxwell–Boltzmann)"
set xlabel "v"
set ylabel "P(v)"
plot \
  "velocity_distributions20.txt" u 1:2 w lp ls 1 t "MC (hist)", \
  "velocity_distributions20.txt" u 1:3 w  l ls 2 t "MB theory"

# Right panel: Position PDF at t = 50
set title "Position PDF at t = 50 (D = 2)"
set xlabel "x"
set ylabel "P(x,t)"
plot \
  "position_distributions20.txt" u 1:2 w lp ls 3 t "MC (hist)", \
  "position_distributions20.txt" u 1:3 w  l ls 4 t "Gaussian  N(0, 2Dt)"

unset multiplot
