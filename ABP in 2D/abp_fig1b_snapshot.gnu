reset
set size ratio -1
set yrange [-1:1]
set xrange [-1:1]
set view map
set title "2D Histogram at t = 1.00" font ",14"
set palette rgbformulae 22,13,-31

set xlabel "x"; set ylabel "y"; set cblabel "P(x,y,t*)"
plot "hist2d_t1.00.dat" matrix with image
