set key left top
set logscale xy
set xlabel "t"
set ylabel "MSD  (\\sigma_x^2, \\sigma_y^2)"


# Load
plot \
  "msd_theory_abp_fig1.txt" u 1:2 w l lw 2 lc "green" t "theory σ_x^2", \
  "msd_theory_abp_fig1.txt" u 1:3 w l lw 2 lc rgb "#1f77b4" t "theory σ_y^2", \
  "msd_sim_abp_fig1.txt"    u 1:2 w p pt 2 ps 0.6 lc "dark-green" t "sim σ_x^2", \
  "msd_sim_abp_fig1.txt"    u 1:3 w p pt 2 ps 0.6 lc "blue" t "sim σ_y^2", \
  2*50*x w l dt 2 lc rgb "black" t "2 D_eff t", \
  0.002*x**3 w l dt 3 lc rgb "gray60" t "t^3 guide", \
  0.00002*x**4 w l dt 3 lc rgb "gray60" t "t^4 guide"
