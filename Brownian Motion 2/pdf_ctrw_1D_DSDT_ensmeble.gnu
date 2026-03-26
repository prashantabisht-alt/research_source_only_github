set grid
set xlabel "x"
set ylabel "P(x,n)"
plot "xPDF_dsdt_p50_nsteps10000.txt" u 1:2 w impulse t "MC pmf", \
     "xPDF_exact_p50_nsteps10000.txt" u 1:2 w l lw 2 t "Exact"
