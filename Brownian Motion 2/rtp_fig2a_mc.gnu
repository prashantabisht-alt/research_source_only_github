# Panel (a)
set xlabel "t"; set ylabel "P(0,t)"; set key top right
set grid back
plot \
  'P0_vs_t_Dp10ptxt'  u 1:2 w l lw 3 lc rgb "#8e44ad" t 'D=0.10 (MC)', \
  'P0_vs_t_Dp175ptxt' u 1:2 w l lw 3 lc rgb "#16a085" t 'D=0.175 (MC)', \
  'P0_vs_t_Dp20ptxt'  u 1:2 w l lw 3 lc rgb "#2980b9" t 'D=0.20 (MC)'
