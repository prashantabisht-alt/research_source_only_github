
set xlabel "lag time τ"
set ylabel "C(τ)"

plot "u_autocorr.txt" u 1:2 w lp lw 2 lc rgb "blue" t "C_est (sim)", \
     "u_autocorr.txt" u 1:3 w l lw 2 lc rgb "red"  t "C_theory"
