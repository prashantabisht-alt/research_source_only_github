# ===== DSDT MSD vs exact =====
# Build tag automatically (edit if you changed p or nsteps)
pcent   = 50          # e.g., p = 0.50 -> 50
nsteps  = 10000
tag     = sprintf("p%02d_nsteps%d", pcent, nsteps)
fmsd    = "msd_dsdt_".tag.".txt"

set term qt size 1100,650
set title sprintf("DSDT MSD (p=%.2f, nsteps=%d)", pcent/100.0, nsteps)
set xlabel "step n"
set ylabel "MSD(n)"
set grid back lw 1
set border lw 1.5
set tics out
set key top left box opaque width -2 spacing 1.2

# Styles
set style line 1 lc rgb "#d62728" pt 5 ps 0.9 lw 1.6   # MC
set style line 2 lc rgb "black"   lw 2.6               # exact

plot \
  fmsd u 1:3 w lp ls 1 t "MSD_MC", \
  fmsd u 1:5 w l  ls 2 t "MSD_exact"
