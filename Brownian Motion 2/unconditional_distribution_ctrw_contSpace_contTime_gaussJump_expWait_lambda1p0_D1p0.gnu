file = "jumps.txt"
D = 1.0
lambda = 1.0
b = sqrt(D/lambda)

bw = 0.05
bin(x) = bw*floor(x/bw) + bw/2.0

# total samples
stats file u 1 name "J" nooutput
N = J_records

set term qt 0
unset logscale y
set title sprintf("Unconditional jump PDF (signed); Laplace b=%.3f", b)
set xlabel "jump"
set ylabel "PDF"
set grid
set style fill solid 0.4
set boxwidth bw absolute

laplace_signed(x) = 1.0/(2.0*b) * exp(-abs(x)/b)

plot file u (bin($1)):(1.0/N/bw) smooth freq w boxes t "Histogram PDF", \
     laplace_signed(x) w l lw 2 t "Laplace theory"

