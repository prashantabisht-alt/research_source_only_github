reset
set term qt enhanced
set pm3d at b
set view map
set palette rgbformulae 22,13,-31
set xlabel "x"
set ylabel "v"
set cblabel "P(x,v)"
set title "Joint distribution P(x,v) at Tmax=50,θ=0.61° from Langevin simulation"
unset key

splot "joint_distribution50.txt" using 1:2:3 with pm3d
