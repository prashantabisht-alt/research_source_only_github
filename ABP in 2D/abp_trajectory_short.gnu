set terminal qt size 500,500
set size ratio -1
set title "Short time trajectory (t << 1/DR)"
set xlabel "x"
set ylabel "y"
plot "traj_short.txt" u 2:3 w lp pt 7 ps 1 lc rgb "red" t "trajectory", \
     "traj_short.txt" u 2:3 every ::0::0 w p pt 7 ps 2 lc rgb "green" t "start"
