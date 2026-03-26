set xlabel "x"; set ylabel "P(x, n = N)"
plot "xPDF_dsdt_p50_nsteps10000.txt"   u 1:2 w p pt 7 ps 0.4 t "MC", \
     "xPDF_exact_p50_nsteps10000.txt"  u 1:2 w l lw 2    t "Exact"

