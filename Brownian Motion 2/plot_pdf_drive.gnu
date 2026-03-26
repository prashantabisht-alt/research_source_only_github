set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Probability Distribution at t = Tmax (Free Overdamped)'
set grid

plot 'pdf_undriven.txt' using 1:2 with lines lw 2 lc rgb 'blue' title 'Undriven', \
     'pdf_driven.txt' using 1:2 with lines lw 2 lc rgb 'red' title 'Driven'
