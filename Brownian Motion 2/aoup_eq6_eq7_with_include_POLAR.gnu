set multiplot layout 1,2 title "AOUP Eq. 6 & Eq. 7 (Empirical vs Theory)" font ",14"

# --- Eq. 6: R_u(t) ---
set title "Eq. 6: Velocity Autocorrelation R_u(t)"
set xlabel "t"
set ylabel "R_u(t)"
set key top right
set grid
plot "eq6_ru.txt" using 1:2 with linespoints title "Empirical" lw 2 lc rgb "blue", \
     "eq6_ru.txt" using 1:3 with lines title "Theory" lw 2 lc rgb "red" dt 2

# --- Eq. 7: M(t) ---
set title "Eq. 7: Memory Kernel M(t)"
set xlabel "t"
set ylabel "M(t)"
set key top right
set grid
plot "eq7_M.txt" using 1:2 with linespoints title "Empirical" lw 2 lc rgb "blue", \
     "eq7_M.txt" using 1:3 with lines title "Theory" lw 2 lc rgb "red" dt 2

unset multiplot
set output
