set xlabel "t"
set ylabel "MSD(t)"
set key top left
set grid

plot "msd_dsdt_p50_nsteps10000.txt" u 1:3 w l lw 2 lc rgb "purple" t "MSD_MC", \
     ""                             u 1:5 w l lw 2 lc rgb "green"  dt 2 t "MSD_exact"
