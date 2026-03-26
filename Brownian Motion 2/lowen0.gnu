#===============================================================
# plot_aoup_multiplot.gp
# Multiplot: C(τ) and M(τ) from AOUP simulation
#===============================================================

set term qt title "AOUP: C(tau) & M(tau)" persist
set multiplot layout 2,1 title "AOUP statistics" font ",14"

# --- 1st plot: C(τ)
set grid
set xlabel "{/Symbol t}"
set ylabel "C({/Symbol t})"
set key outside top right
plot "u_autocorr.txt" using 1:2 with points pt 7 ps 0.5 title "C_est", \
     ""                using 1:3 with lines  lw 2    title "C_theory"

# --- 2nd plot: M(τ)
set grid
set xlabel "{/Symbol t}"
set ylabel "M({/Symbol t})"
set key outside top right
plot "M_kernel.txt" using 1:2 with points pt 7 ps 0.5 title "M_est", \
     ""              using 1:3 with lines  lw 2    title "M_theory"

unset multiplot
