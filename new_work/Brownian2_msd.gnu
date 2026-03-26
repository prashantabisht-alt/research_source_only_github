set term qt size 900,550
set xlabel "time t"
set ylabel "<x^2(t)>"
set grid
set key top left

D = 1.0   # same as in your Fortran code

plot \
  "Brownian_msd.txt" u 1:2 w lp pt 5 ps 1.0 lw 1.5 lc rgb "red"   t "MSD (ensemble)", \
  2*D*x              w l  dt 2 lw 2   lc rgb "black" t "theory: 2Dt"
