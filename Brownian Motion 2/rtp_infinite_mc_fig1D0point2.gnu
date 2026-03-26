set xlabel "x"
set ylabel "P(x,t) for D=0.2"
set key left top

times = "tp2 tp4 tp6 tp8 t1p0 t1p2 t1p4 t1p6 t1p8 t2p0"
plot for [ti in times] sprintf("P_Dp20_%sptxt", ti) u 1:2 w lp lw 2 t ti
