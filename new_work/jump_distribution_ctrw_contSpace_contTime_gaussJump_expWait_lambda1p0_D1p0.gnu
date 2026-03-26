file = "jumps.txt"   # cols: jump  dt
D    = 1.0

# standard normal pdf
pi = 3.141592653589793
phi(x) = 1.0/sqrt(2*pi) * exp(-0.5*x*x)

# binning for standardized z
bw = 0.1
bin(x) = bw*floor(x/bw) + bw/2.0

# count records
stats file u 1 name "J" nooutput
N = J_records

set term qt 0
set title "Standardized jumps z = jump / sqrt(2 D Δt)  vs  N(0,1)"
set xlabel "z"
set ylabel "PDF"
set grid
set style fill solid 0.4
set boxwidth bw absolute

plot file u ( bin($1/sqrt(2*D*$2)) ):(1.0/N/bw) smooth freq w boxes t "Histogram of z", \
     phi(x) w l lw 2 t "N(0,1)"
