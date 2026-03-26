set term qt size 1100,650
set title "CS-CTRW MSD(t) vs 2Dt"
set xlabel "t"; set ylabel "MSD(t)"
set grid back lw 1; set border lw 1.5; set tics out
set key top left box opaque
D=1.0
plot "msd_ensemble_csctrw_lambda1p0_D1p0_tmax1e4.txt" u 1:2 w lp pt 5 ps 0.9 lw 1.6 lc rgb "#d62728" t "MSD (ensemble)", \
     2*D*x w l lw 2.6 lc rgb "black" t "Theory 2Dt"
