set grid
set xlabel 'time t'
set xrange [0:50]
set ylabel 'site (wrapped: 0..L-1)'
plot 'ctrw_single_particle.txt' using 1:2 with steps lw 1.5 title 'wrapped pos'
