reset
set term qt enhanced
set view map
set pm3d at b
set palette rgbformulae 22,13,-31
set xlabel "x"
set ylabel "v"
set cblabel "P(x,v)"
unset key
set multiplot layout 2,2 title "Joint distribution theory  P(x,v) — pm3d heatmaps"

# Panel 1: t = 5
set title "T = 5"
splot "joint_theory5.txt" using 1:2:3 with pm3d

# Panel 2: t = 20
set title "T = 20"
splot "joint_theory20.txt" using 1:2:3 with pm3d

# Panel 3: t = 50
set title "T = 50"
splot "joint_theory50.txt" using 1:2:3 with pm3d

# Panel 4: t = 100
set title "T = 100"
splot "joint_theory100.txt" using 1:2:3 with pm3d

unset multiplot
