reset
set term qt size 1200,1200
set multiplot layout 2,2 title "Joint P(x,v) contours — Simulation vs Theory"

# ===== Common setup =====
set view map
unset surface
unset clabel
set contour base
set grid front lw 0.5
set key top right box opaque
set xrange [-50:50]
set yrange [-5:5]

# Contour levels
set cntrparam levels discrete 0.002,0.004,0.006,0.008,0.010,0.012

# Function to add a contour label
ContourLabel() = "Contour levels: 0.002  0.004  0.006  0.008  0.010  0.012"

# --- Panel 1: Tmax = 5 ---
set title "Tmax = 5 (θ = 11.06°)"
set xlabel "x"; set ylabel "v"
set label 1 ContourLabel() at graph 0.02, graph 0.92 font "Helvetica,10"
set table "sim_contours5.dat"
splot "joint_distribution5.txt" u 1:2:3
unset table
set table "theory_contours5.dat"
splot "joint_theory5.txt" u 1:2:3
unset table
plot "sim_contours5.dat" u 1:2 w l lc rgb "black" lw 1.2 t "Simulation", \
     "theory_contours5.dat" u 1:2 w l lc rgb "red" lw 1.6 t "Theory"
unset label 1

# --- Panel 2: Tmax = 20 ---
set title "Tmax = 20 (θ = 1.77°)"
set xlabel "x"; set ylabel "v"
set label 1 ContourLabel() at graph 0.02, graph 0.92 font "Helvetica,10"
set table "sim_contours20.dat"
splot "joint_distribution20.txt" u 1:2:3
unset table
set table "theory_contours20.dat"
splot "joint_theory20.txt" u 1:2:3
unset table
plot "sim_contours20.dat" u 1:2 w l lc rgb "black" lw 1.2 t "Simulation", \
     "theory_contours20.dat" u 1:2 w l lc rgb "red" lw 1.6 t "Theory"
unset label 1

# --- Panel 3: Tmax = 50 ---
set title "Tmax = 50 (θ = 0.61°)"
set xlabel "x"; set ylabel "v"
set label 1 ContourLabel() at graph 0.02, graph 0.92 font "Helvetica,10"
set table "sim_contours50.dat"
splot "joint_distribution50.txt" u 1:2:3
unset table
set table "theory_contours50.dat"
splot "joint_theory50.txt" u 1:2:3
unset table
plot "sim_contours50.dat" u 1:2 w l lc rgb "black" lw 1.2 t "Simulation", \
     "theory_contours50.dat" u 1:2 w l lc rgb "red" lw 1.6 t "Theory"
unset label 1

# --- Panel 4: Tmax = 100 ---
set title "Tmax = 100 (θ = 0.29°)"
set xlabel "x"; set ylabel "v"
set label 1 ContourLabel() at graph 0.02, graph 0.92 font "Helvetica,10"
set table "sim_contours100.dat"
splot "joint_distribution100.txt" u 1:2:3
unset table
set table "theory_contours100.dat"
splot "joint_theory100.txt" u 1:2:3
unset table
plot "sim_contours100.dat" u 1:2 w l lc rgb "black" lw 1.2 t "Simulation", \
     "theory_contours100.dat" u 1:2 w l lc rgb "red" lw 1.6 t "Theory"
unset label 1

unset multiplot
