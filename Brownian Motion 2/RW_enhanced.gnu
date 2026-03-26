dtmin=0
dtmax=10
nbins=200
dx=(dtmax-dtmin)/nbins
stats "waiting_time_hist_counts.txt using ($1*$2) name "M" nooutput
mean_dt = M_sum*dx
lamda_hat=1.0/mean_dt
print sprintf("λ : %.6f", lamda_hat)
f_pdf(x) = lamda_hat*exp(-lamda_hat*x)
set title "delta t with fitted exponentIAL)
set xlabel "delta t"
set ylabel "PDF"
plot "waiting_time_hist_counts.txt u 1:2 w lp lw 2 title "simulation_PDF", /f_pdf(x) with lines lw 2 
