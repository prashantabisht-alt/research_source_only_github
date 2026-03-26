set term qt size 1000,600
set xlabel "time t"
set ylabel "position x(t)"
set grid
set key outside right spacing 1.2

# One substream per trajectory: only rows with $2==k (traj id), print t and x
plot for [k=1:10] sprintf("< awk '$2==%d {print $1, $3}' Brownian_trajectories.txt", k) \
     u 1:2 w lp pt (4+2*k) ps 0.7 lw 1.6 t sprintf("traj %d",k)
