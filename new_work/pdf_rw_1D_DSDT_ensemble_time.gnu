set xlabel "x"
set ylabel "P(x,t=1000)"
set key top right
set grid
plot "xPDF_dsdt_p50_nsteps10000_t1000p000.txt" u 1:2 w p pt 7 ps 0.4 lc rgb "purple" t "MC", \
     "xPDF_exact_p50_nsteps10000_t1000p000.txt" u 1:2 w p pt 7 ps 0.6 lc rgb "green"  t "Exact"

