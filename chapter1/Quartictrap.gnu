# ===== Quartic trap: P(x) with Gaussian noise =====
set term qt size 1100,650
set title "Quartic trap: Position PDF (Gaussian noise)"
set xlabel "x"
set ylabel "P(x)"
set grid back lw 1
set border lw 1.6
set tics out
set key top right box opaque width -2 spacing 1.2
set samples 800

# Parameters
a   = 0.1
kBT = 1.0
pi  = 3.141592653589793

# Normalization constant (Gamma(1/4))
Z = 0.5*gamma(0.25)*(kBT/a)**0.25
p_quart(x) = exp(-a*x**4/kBT)/Z

# Styles
set style line 3 lc rgb "#1f77b4" pt 7 ps 0.9 lw 1.6   # MC: blue triangles
set style line 4 lc rgb "black"   lw 2.8               # Theory

# --- Main plot ---
set multiplot

plot \
  "pdf_quart_gauss_full.txt" u 1:2 w lp ls 3 t "MC: Gaussian noise", \
  p_quart(x) w l ls 4 t "Boltzmann: ∝ exp(-a x^4/kBT)"

# --- Inset plot (log-scale y) ---
set origin 0.55,0.45     # inset position (tune)
set size   0.4,0.4       # inset size (tune)
set logscale y
unset key
unset title
set xlabel ""
set ylabel ""
set grid back lw 0.8
set tics out

replot

unset multiplot

