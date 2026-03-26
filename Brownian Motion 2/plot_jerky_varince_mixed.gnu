set terminal pngcairo size 1400,1000 enhanced font 'Arial,14'
set key top left
GAUSS = "mixed_stats.txt"
EXPN  = "mixed_stats_exponential.txt"
set multiplot layout 1,3 title "Variance vs Time (log y) — Gaussian vs Exponential" font ",16"

set xlabel "t"; set ylabel "Var(a)"
plot GAUSS u 1:3 w l lw 2 lc rgb "blue" title "Gaussian", \
     EXPN  u 1:3 w l lw 2 lc rgb "red"  title "Exponential"

set xlabel "t"; set ylabel "Var(v)"
plot GAUSS u 1:5 w l lw 2 lc rgb "blue" title "Gaussian", \
     EXPN  u 1:5 w l lw 2 lc rgb "red"  title "Exponential"

set xlabel "t"; set ylabel "Var(x)"
plot GAUSS u 1:7 w l lw 2 lc rgb "blue" title "Gaussian", \
     EXPN  u 1:7 w l lw 2 lc rgb "red"  title "Exponential"

unset multiplot
replot
set output
