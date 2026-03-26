set term qt size 1100,650
set title "MSD(t): CTRW (exp waits, exact) vs CSRW (discrete-time)"
set xlabel "t"
set ylabel "MSD(t)"
set grid back lw 1
set border lw 1.5
set tics out
set key top left box opaque width -2 spacing 1.2

D = 1.0

# styles
set style line 1 lc rgb "#d62728" pt 5 ps 0.9 lw 1.6  # CTRW = red squares + line
set style line 2 lc rgb "#1f77b4" pt 7 ps 0.9 lw 1.6  # CSRW = blue triangles + line
set style line 3 lc rgb "black"   lw 2.6               # theory = solid black

plot \
  "msd_ensemble.txt" u 1:2 w lp ls 1 t "CTRW (exp waits)", \
  "msd_ensemble_csrw.txt" u 1:2 w lp ls 2 t "CSRW (discrete-time)", \
  2*D*x                   w  l  ls 3 t "Theory: 2Dt"
