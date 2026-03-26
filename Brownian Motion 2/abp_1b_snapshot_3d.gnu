reset
set term qt
set xlabel "x"; set ylabel "y"; set zlabel "P(x,y)"
set hidden3d; set pm3d depthorder
splot "hist2d_t1.00.dat" matrix w pm3d notitle
