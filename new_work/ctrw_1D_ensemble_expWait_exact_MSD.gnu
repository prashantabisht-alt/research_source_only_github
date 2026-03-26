set title "Ensemble MSD vs t"
set xlabel "t"
set ylabel "MSD(t)"
set grid
D = 1.0
plot "msd_ensemble.txt" u 1:2 w lp lw 2.5 t "Ensemble MSD", \
     2*D*x w l dt 2 lw 2 t "2 D t"
