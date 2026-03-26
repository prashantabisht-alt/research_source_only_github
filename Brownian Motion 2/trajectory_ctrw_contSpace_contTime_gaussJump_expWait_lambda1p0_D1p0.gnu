set xlabel "Time"
set ylabel "Position"
plot "< head -n 500 trajectory_ctrw_contSpace_contTime_gaussJump_expWait_lambda1p0_D1p0.txt" using 1:2 with lp lw 2  title "First 500 jumps"

