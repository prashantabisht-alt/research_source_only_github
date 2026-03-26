reset
set term qt enhanced
set view map
set pm3d at b
set palette rgbformulae 22,13,-31

set xlabel "x"
set ylabel "v"
set cblabel "P(x,v)"
set cbrange [0:0.012]
set multiplot layout 2,2 title "Joint distribution P(x,v) — pm3d heatmaps"

# Panel 1: t = 50
set title "T = 50"
splot "joint_distribution50.txt" using 1:2:3 with pm3d
# Panel 1: t = 50
set title "T = 50"
splot "joint_theory50.txt" using 1:2:3 with pm3d

unset multiplot
