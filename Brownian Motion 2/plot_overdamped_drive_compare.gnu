set multiplot layout 2,1 title "Overdamped Langevin: Undriven vs Driven"

# Panel 1: Mean Position
set xlabel 'Time t (s)'
set ylabel '<x(t)>'
set title 'Mean Position Response'
set grid
plot 'mean_pos_compare.txt' using 1:2 with lines lw 2 lc rgb 'blue' title 'Undriven', \
     'mean_pos_compare.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Driven'

# Panel 2: MSD
set xlabel 'Time t (s)'
set ylabel 'MSD <x^2(t)>'
set title 'Mean Squared Displacement'
set grid
plot 'msd_compare.txt' using 1:2 with lines lw 2 lc rgb 'blue' title 'Undriven', \
     'msd_compare.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Driven'

unset multiplot
