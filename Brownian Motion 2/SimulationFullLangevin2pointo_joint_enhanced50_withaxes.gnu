reset
set term qt
set view map
unset surface; unset clabel
set contour base
set cntrparam levels discrete 0.002,0.004,0.006,0.010

set table "sim_contours.dat";    splot "joint_distribution50.txt" u 1:2:3; unset table
set table "theory_contours.dat"; splot "joint_theory50.txt"       u 1:2:3; unset table

plot "sim_contours.dat"    u 1:2 w l lc rgb "black" lw 1.2 t "sim contours", \
     "theory_contours.dat" u 1:2 w l lc rgb "red"   lw 1.6 t "theory contours", \
     "axes_theta_theory50.txt" u 1:2 w l lc rgb "orange" lw 2 t "theory axes", \
     "axes_theta_sim50.txt"    u 1:2 w l lc rgb "blue"   dt 2 lw 1.6 t "sim axes", \
     0 w l lc rgb "gray" dt 2 lw 1.0 t "x=0", \
     "" u ($1):(0) w l lc rgb "gray" dt 2 lw 1.0 t "v=0"
