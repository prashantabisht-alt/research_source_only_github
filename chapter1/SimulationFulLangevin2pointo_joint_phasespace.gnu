reset
set term qt enhanced
set xlabel "x"
set ylabel "v"
set title "Phase Space at t = T_max"
set grid
plot "phase_space_points.txt" using 1:2 with dots notitle
