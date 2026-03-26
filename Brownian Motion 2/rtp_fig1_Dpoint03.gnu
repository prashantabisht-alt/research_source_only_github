reset
set title "RTP Fig.1 — D=0.03"
set xlabel "x"
set ylabel "P(x,t)"
set grid
set logscale y
set key top left

# columns: 1=x, 2..11 = t=0.2,0.4,...,2.0
plot for [c=2:11] "fig1_D0.03.txt" u 1:c w l lw 2 t sprintf("t=%.1f", 0.2*(c-1))
