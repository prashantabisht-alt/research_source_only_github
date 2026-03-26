reset
set term qt size 900,400

# max from simulation panel
set cbrange [0:0.014]

set multiplot layout 1,2 title "Joint P(x,v) heatmaps  T=50  (θ = 0.61°)"

# --- Common setup ---
set palette rgbformulae 22,13,-31
set xlabel "x"
set ylabel "v"
set cblabel "P(x,v)"
set xrange [-50:50]
set yrange [-5:5]

# --- Simulation heatmap ---
set title "Simulation"
plot "joint_distribution50.txt" u 1:2:3 with image notitle

# --- Theory heatmap ---
set title "Theory"
plot "joint_theory50.txt" u 1:2:3 with image notitle

unset multiplot
