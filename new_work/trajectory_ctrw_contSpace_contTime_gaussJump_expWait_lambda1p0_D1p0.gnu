set term qt size 1000,600
set title "CTRW trajectory – first 50 jumps"
set xlabel "Time"
set ylabel "Position x(t)"
set grid back lw 1
set border lw 1.5
set tics out
set key top left box opaque

# Style: red squares with connecting line
set style line 1 lc rgb "#d62728" pt 5 ps 1.0 lw 2

plot "< head -n 50 trajectory_ctrw_contSpace_contTime_gaussJump_expWait_lambda1p0_D1p0.txt" \
     u 1:2 w lp ls 1 t "First 50 jumps"

