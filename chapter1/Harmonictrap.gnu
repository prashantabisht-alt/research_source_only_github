set term qt size 1100,650
set title "Harmonic trap: Position PDF (Gaussian noise)"
set xlabel "x"
set ylabel "P(x)"
set grid back lw 1
set border lw 1.6
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 800


k   = 0.2
kBT = 1.0
sigma2 = kBT/k
pi = 3.141592653589793
p_harm(x) = 1.0/sqrt(2*pi*sigma2) * exp(-x*x/(2*sigma2))

# Styles
set style line 1 lc rgb "#d62728" pt 5 ps 0.9 lw 1.6   # MC: red squares
set style line 2 lc rgb "black"   lw 2.8               # Theory


set multiplot

plot \
  "pdf_harm_gauss_full.txt" u 1:2 w lp ls 1 t "MC: Gaussian noise", \
  p_harm(x) w l ls 2 t "Boltzmann: N(0, kBT/k)"

set origin 0.60,0.45     # lower-left corner of inset (in relative coords 0–1)
set size   0.4,0.4       # width, height of inset
set logscale y
unset key
set title ""
set xlabel ""
set ylabel ""
set tics out
set grid back lw 0.8

replot

unset multiplot

