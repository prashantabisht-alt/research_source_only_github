set xlabel "n"; set ylabel "MSD"
plot "msd_dsdt_p50_nsteps10000.txt" u 1:3 w l lw 2 t "MSD_MC", \
     ""                              u 1:5 w l lw 2 dt 2 t "MSD_exact"
