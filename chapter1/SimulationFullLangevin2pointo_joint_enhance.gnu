reset
set term qt enhanced
set view map            # 2D color map (top view)
set pm3d at b           # pm3d coloring at bottom surface
set palette rgbformulae 22,13,-31  # nice rainbow-like palette
set xlabel "x"
set ylabel "v"
set cblabel "P(x,v)"
set title "Joint distribution P(x,v) at Tmax=20 ,θ=1.77° from Langevin simulation"
unset key
splot "joint_distribution20.txt" using 1:2:3 with pm3d
