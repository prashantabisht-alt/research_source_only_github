set xlabel "lag time τ"
set ylabel "M(τ)"

plot "M_kernel.txt" u 1:2 w lp pt 2 ps 1.2 lc rgb "blue"  t "M(simulation)", \
     "M_kernel.txt" u 1:3 w l  lw 2   lc rgb "dark-red" t "(M_theory)"
