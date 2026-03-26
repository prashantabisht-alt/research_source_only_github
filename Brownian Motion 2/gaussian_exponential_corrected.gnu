set title "PDF of Final Positions (Log Scale)"
set xlabel "x"
set ylabel "log P(x)"
set logscale y
plot \
    "gaussian_simulation1.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "blue" title "Gaussian Sim", \
    "exponential_corrected1.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "red" title "Exponential Sim", \
    (1.0/sqrt(4.0*pi*1.0*10000*0.001)) * exp(-x**2/(4.0*1.0*10000*0.001)) lc rgb "black" lw 2 title "Theory (Gaussian)"
