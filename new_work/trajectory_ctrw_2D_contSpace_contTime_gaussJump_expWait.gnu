file = "trajectory_ctrw_2D_contSpace_contTime_gaussJump_expWait_lambda1p0_D1p0.txt"
set term qt
set title "2D CTRW Trajectory (colored by time)"
set xlabel "x"
set ylabel "y"
set grid
set palette rgb 33,13,10   # nice gradient
plot file using 2:3:1 with lines lc palette title "Path"
