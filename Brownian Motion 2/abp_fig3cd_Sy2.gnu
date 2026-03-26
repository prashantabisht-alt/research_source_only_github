set logscale x
set logscale y
set xlabel "u = t D_R"
set ylabel "(v_0/(y_0 D_R))^{1/6} S_y"
set key bottom left
set grid y

# slope guides (paper-style)
RAPu(x)  = x**(-0.25)   # u^{-1/4}
DIFFu(x) = x**(-0.50)   # u^{-1/2}

plot \
  "Sy_collapsed_DR0.0005.txt" u 1:2 w l lw 2 t "D_R = 0.0005", \
  "Sy_collapsed_DR0.0010.txt" u 1:2 w l lw 2 t "D_R = 0.0010", \
  "Sy_collapsed_DR0.0050.txt" u 1:2 w l lw 2 t "D_R = 0.0050", \
  RAPu(x)  w l dt 3 lw 1.6 t "u^{-1/4} (RAP)", \
