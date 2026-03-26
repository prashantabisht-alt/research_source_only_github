set term qt enhanced
set multiplot layout 1,3 title "RTP in box (ell=1): steady state vs analytic"
set xlabel "x"; set ylabel "P(x)"; set xrange [-1:1]
# set logscale y      # uncomment to match paper's wall emphasis

plot 'fig3_D0p10ptxt' u 1:2 w p pt 7 ps 0.25 lc rgb "black" t 'MC', \
     'fig3_D0p10ptxt' u 1:3 w l lw 3 lc rgb "orange" t 'analytic'
plot 'fig3_D0p50ptxt' u 1:2 w p pt 7 ps 0.25 lc rgb "black" t 'MC', \
     'fig3_D0p50ptxt' u 1:3 w l lw 3 lc rgb "orange" t 'analytic'
plot 'fig3_D1p00ptxt' u 1:2 w p pt 7 ps 0.25 lc rgb "black" t 'MC', \
     'fig3_D1p00ptxt' u 1:3 w l lw 3 lc rgb "orange" t 'analytic'
unset multiplot
