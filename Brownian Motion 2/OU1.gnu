set title "PDF: Histogram vs Gaussian"
set xlabel "u"
set ylabel "Probability density"

plot \
    "ou_hist_vs_gauss.txt" using 1:2 with lines lw 2 title "P_sim(u)", \
    "" using 1:3 with lines lw 2 dt 2 title "Gaussian theory"
