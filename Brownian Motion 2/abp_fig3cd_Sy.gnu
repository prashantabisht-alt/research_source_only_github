set logscale x
set logscale y
set xlabel "t"
set ylabel "S_y(t; y_0)"
set key bottom left

RAP(x)  = x**(-0.25)
DIFF(x) = x**(-0.50)

plot \
 "Sy_DR0.0005.txt" u 1:2 w l lw 2 t "D_R=0.0005", \
 "Sy_DR0.0010.txt" u 1:2 w l lw 2 t "D_R=0.0010", \
 "Sy_DR0.0050.txt" u 1:2 w l lw 2 t "D_R=0.0050", \
 "Sy_RAP_asym_DR0.0005.txt" u 1:2 w l dt 2 lw 2 t "RAP asymptote (0.0005)", \
 "Sy_RAP_asym_DR0.0010.txt" u 1:2 w l dt 2 lw 2 t "RAP asymptote (0.0010)", \
 "Sy_RAP_asym_DR0.0050.txt" u 1:2 w l dt 2 lw 2 t "RAP asymptote (0.0050)", \
 RAP(x)  w l dt 3 lw 1.6 t "t^{-1/4}", \
 DIFF(x) w l dt 4 lw 1.6 t "t^{-1/2}"
