set xlabel 'Position x'
set ylabel 'P(x)'
set title 'Probability Distribution with OU Noise (Overdamped)'
set grid

plot 'pdf_ou.txt' using 1:2 with lines lw 2 title sprintf("OU Noise, tau_c=%.2f", 0.1)
