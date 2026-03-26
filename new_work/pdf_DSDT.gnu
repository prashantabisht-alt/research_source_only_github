# ===== DSDT PMF at final time: MC vs Exact =====
pcent   = 50
nsteps  = 10000
tag     = sprintf("p%02d_nsteps%d", pcent, nsteps)
fmc     = "xPDF_dsdt_".tag.".txt"
fex     = "xPDF_exact_".tag.".txt"

set term qt size 1100,650
set title sprintf("DSDT P(X_n) at n=%d (p=%.2f)", nsteps, pcent/100.0)
set xlabel "x (lattice site)"
set ylabel "P(X_n = x)"
set grid back lw 1
set border lw 1.8
set tics out
set key top right box opaque width -2 spacing 1.2

# Styles
set style line 11 lc rgb "#1f77b4" pt 5 ps 0.9 lw 1.6   # MC: blue squares + line
set style line 12 lc rgb "#d62728" lw 2.4              # Exact: dark red stems

plot \
  fmc u 1:2 w lp       ls 11 t "MC (simulation)", \
  fex u 1:2 w impulses ls 12 t "Exact (binomial)"
