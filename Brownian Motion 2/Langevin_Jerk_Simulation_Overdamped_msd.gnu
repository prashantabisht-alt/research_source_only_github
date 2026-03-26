set title 'MSD vs. Time (Third-Order Lambda)'
set xlabel 'Time t (s)'
set ylabel 'MSD <x^2> (units^2)'
f_msd(x) = 2.0 * x  #Approximate long-time behavior
plot 'msd_third_order_lambda.txt' using 1:2 with points title 'Simulated', \
     f_msd(x) with lines title 'Theoretical 2Dt'
