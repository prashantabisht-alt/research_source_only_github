set xlabel "lag time τ"
set ylabel "M(τ)"

plot "M_kernel.txt" u 1:2 w l lw 2 lc rgb "blue" t "M_est (sim)", \
     "M_kernel.txt" u 1:3 w l lw 2 lc rgb "red"  t "M_theory"
