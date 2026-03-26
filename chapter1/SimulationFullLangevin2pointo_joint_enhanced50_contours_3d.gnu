reset
set term qt
set hidden3d
set contour base
set clabel
set cntrparam levels discrete 0.002,0.004,0.006,0.008,0.010,0.012

set xlabel "x"
set ylabel "v"
set zlabel "P(x,v)"

# Simulation contours as surface
splot "joint_distribution50.txt" u 1:2:3 w l lc rgb "black" lw 1.2 t "Simulation", \
      "joint_theory50.txt"       u 1:2:3 w l lc rgb "red"   lw 1.6 t "Theory"
