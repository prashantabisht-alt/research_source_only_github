times = "tp2 tp4 tp6 tp8 t1p0 t1p2 t1p4 t1p6 t1p8 t2p0"
plot \
  for [ti in times] sprintf("P_eq10_Dp03_%sptxt", ti) u 1:2 w l lw 2 t sprintf("analytic %s", ti), \
  for [ti in times] sprintf("P_Dp03_%sptxt", ti)      u 1:2 w p pt 6 ps 0.3 t sprintf("MC %s", ti)
