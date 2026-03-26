set xlabel 'Time t (s)'
set ylabel '<v(t)v(0)>'
set title 'Velocity Autocorrelation: Undriven vs Driven (Overdamped)'
set grid

plot 'vacf_compare.txt' using 1:2 with lines lw 2 lc rgb 'blue' title 'Undriven', \
     'vacf_compare.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Driven'
