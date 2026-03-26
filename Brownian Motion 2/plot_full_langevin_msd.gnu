set xlabel 'x'
set ylabel 'MSD <x^2> (units^2)'
set title 'MSD for Full Langevin (x = 50)' font ",14"
set key bottom right
set grid

# --- Theoretical Curves ---
f_free(x)  = 2.0 * x
f_harm(x)  = 5.0 * (1.0 - exp(-x / 2.0))

# --- Plot ---
plot \
  'msd_free_gauss.txt' using 1:2 with points pt 7 ps 0.6 lc 'blue'    title 'Free (Gaussian)', \
  'msd_free_exp.txt'   using 1:2 with points pt 7 ps 0.6 lc 'red'     title 'Free (Exponential)', \
  f_free(x)            with lines      lw 2  lc 'black'               title 'Theoretical Free', \
  'msd_harm_gauss.txt' using 1:2 with points pt 6 ps 0.6 lc 'cyan'    title 'Harmonic (Gaussian)', \
  'msd_harm_exp.txt'   using 1:2 with points pt 6 ps 0.6 lc 'magenta' title 'Harmonic (Exponential)', \
  f_harm(x)            with lines      lw 2  lc 'gray'                title 'Theoretical Harmonic'
