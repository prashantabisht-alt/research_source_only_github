reset

set view map
set size ratio -1
set xlabel "x"; set ylabel "y"; set cblabel "P(x,y,t*)"
plot "hist2d_t1.00.dat" matrix with image
