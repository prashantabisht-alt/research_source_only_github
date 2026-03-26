reset
set title "RTP Fig.1 — D=0.20"
set xlabel "x"
set ylabel "P(x,t)"
set grid
set logscale y
set key top left

plot for [c=2:11] "fig1_D0.20.txt" u 1:c w l lw 2 t sprintf("t=%.1f", 0.2*(c-1))

