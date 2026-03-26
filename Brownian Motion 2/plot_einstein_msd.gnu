set xlabel 'Time t (s)'
set xrange [0:0.2]
set ylabel 'MSD <x^2> (units^2)'
set title 'Einstein Relation Verification (Overdamped, t = 10)'
set key bottom right
set grid

f(x) = 2.0 * x

plot \
  'msd_gauss_einstein.txt' using 1:2 with points pt 2 ps 0.8  lc 'blue' title 'Gaussian Noise', \
  'msd_exp_einstein.txt' using 1:2 with points pt 2 ps 0.8  lc 'red' title 'Exponential Noise', \
   f(x) with lines lw 2 lc 'black' title 'Theoretical 2Dt'
