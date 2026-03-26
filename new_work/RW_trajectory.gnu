set term qt size 1000,600
set grid
set xlabel 'time t'
set ylabel 'site (wrapped: 0..L-1)'
set xrange [0:50]
set yrange [0:999]   # since L = 1000

plot "< awk '$1<=50 {print $1, $2}' ctrw_single_particle.txt" \
     u 1:2 w steps lw 1.8 lc rgb "blue" t 'wrapped position'
