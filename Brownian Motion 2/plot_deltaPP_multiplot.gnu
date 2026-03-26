set multiplot layout 2,2 title "δ''-forcing (Route A): PDFs and moments"

# Acceleration PDF
set title "Acceleration PDF (simulation)"
set xlabel "a"
set ylabel "P(a)"
plot "acceleration_dist_pp.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "navy" title "sim"

# Velocity PDF
set title "Velocity PDF (simulation)"
set xlabel "v"
set ylabel "P(v)"
plot "velocity_dist_pp.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "dark-green" title "sim"

# Position PDF
set title "Position PDF (simulation)"
set xlabel "x"
set ylabel "P(x)"
plot "position_dist_pp.txt" using 1:2 with points pt 7 ps 0.5 lc rgb "dark-red" title "sim"

# Moments vs time
set title "Running means/variances"
set xlabel "t"
set ylabel "value"
plot \
  "jerkpp_stats.txt" using 1:2 with lines lw 2 lc rgb "navy" title "<a>", \
  "" using 1:3 with lines lw 2 lc rgb "skyblue" title "Var(a)", \
  "" using 1:4 with lines lw 2 lc rgb "dark-green" title "<v>", \
  "" using 1:5 with lines lw 2 lc rgb "green" title "Var(v)", \
  "" using 1:6 with lines lw 2 lc rgb "dark-red" title "<x>", \
  "" using 1:7 with lines lw 2 lc rgb "red" title "Var(x)"

unset multiplot
set output

