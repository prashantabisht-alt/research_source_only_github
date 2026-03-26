set terminal qt size 500,500
set size ratio -1
set title "Long time trajectory (t >> 1/DR)"
set xlabel "x"
set ylabel "y"
plot "traj_long.txt" u 2:3 w lp pt 7 ps 0.6 lc rgb "red" t "trajectory", \
     "traj_long.txt" u 2:3 every ::0::0 w p pt 7 ps 2 lc rgb "green" t "start"
