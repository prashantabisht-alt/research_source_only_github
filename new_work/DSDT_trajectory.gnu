# ===== DSDT: first 50 steps with markers =====
set term qt size 1000,600
set title "DSDT: First 50 steps of one trajectory"
set xlabel "Step n"
set ylabel "X_n"
set grid back lw 1
set border lw 1.5
set tics out
set key off

# Styles
set style line 1 lc rgb "#1f77b4" lw 2       # staircase line
set style line 2 lc rgb "#d62728" pt 7 ps 1.2 # red circles for points

# Plot only first 51 rows (n=0…50)
plot "< head -n 51 traj1.txt" u 1:2 w steps ls 1, \
     "< head -n 51 traj1.txt" u 1:2 w p ls 2
