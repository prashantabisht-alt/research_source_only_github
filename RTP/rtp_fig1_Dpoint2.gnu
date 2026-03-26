set term qt size 1100,650
set title "RTP on infinite line — D = 0.20"
set xlabel "x"
set ylabel "P(x,t)"
set xrange [-2.2:2.2]
set yrange [0:*]
set grid ytics
set key top right spacing 1.15

set style line 101 lc rgb "#1f77b4" pt 7  ps 0.8 lw 1
set style line 201 lc rgb "#1f77b4" lw 2
set style line 102 lc rgb "#ff7f0e" pt 5  ps 0.9 lw 1
set style line 202 lc rgb "#ff7f0e" lw 2
set style line 103 lc rgb "#2ca02c" pt 9  ps 0.8 lw 1
set style line 203 lc rgb "#2ca02c" lw 2
set style line 104 lc rgb "#d62728" pt 13 ps 0.9 lw 1
set style line 204 lc rgb "#d62728" lw 2
set style line 105 lc rgb "#9467bd" pt 11 ps 0.9 lw 1
set style line 205 lc rgb "#9467bd" lw 2
set style line 106 lc rgb "#8c564b" pt 4  ps 0.9 lw 1
set style line 206 lc rgb "#8c564b" lw 2
set style line 107 lc rgb "#e377c2" pt 6  ps 0.9 lw 1
set style line 207 lc rgb "#e377c2" lw 2
set style line 108 lc rgb "#7f7f7f" pt 8  ps 0.9 lw 1
set style line 208 lc rgb "#7f7f7f" lw 2
set style line 109 lc rgb "#bcbd22" pt 12 ps 0.9 lw 1
set style line 209 lc rgb "#bcbd22" lw 2
set style line 110 lc rgb "#17becf" pt 10 ps 1.0 lw 1
set style line 210 lc rgb "#17becf" lw 2

plot \
  "fig1_D0.20.txt"        u 1:2  w p ls 101 t "t=0.2 (sim)", \
  ""                      u 1:3  w p ls 102 t "t=0.4 (sim)", \
  ""                      u 1:4  w p ls 103 t "t=0.6 (sim)", \
  ""                      u 1:5  w p ls 104 t "t=0.8 (sim)", \
  ""                      u 1:6  w p ls 105 t "t=1.0 (sim)", \
  ""                      u 1:7  w p ls 106 t "t=1.2 (sim)", \
  ""                      u 1:8  w p ls 107 t "t=1.4 (sim)", \
  ""                      u 1:9  w p ls 108 t "t=1.6 (sim)", \
  ""                      u 1:10 w p ls 109 t "t=1.8 (sim)", \
  ""                      u 1:11 w p ls 110 t "t=2.0 (sim)", \
  "fig1_theory_D0.20.txt" u 1:2  w l ls 201 t "t=0.2 (theory)", \
  ""                      u 1:3  w l ls 202 notitle, \
  ""                      u 1:4  w l ls 203 notitle, \
  ""                      u 1:5  w l ls 204 notitle, \
  ""                      u 1:6  w l ls 205 notitle, \
  ""                      u 1:7  w l ls 206 notitle, \
  ""                      u 1:8  w l ls 207 notitle, \
  ""                      u 1:9  w l ls 208 notitle, \
  ""                      u 1:10 w l ls 209 notitle, \
  ""                      u 1:11 w l ls 210 notitle


