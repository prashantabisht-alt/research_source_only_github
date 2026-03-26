
dtmin = 0.0
dtmax = 10.0
nbins = 200
dx    = (dtmax - dtmin) / nbins


stats "waiting_time_hist_pdf.txt" using ($1*$2) name "M" nooutput
mean_dt   = M_sum * dx
lambda_hat = 1.0 / mean_dt
print sprintf("Estimated λ̂ from histogram PDF: %.6f", lambda_hat)


f_pdf(x) = lambda_hat * exp(-lambda_hat * x)


set term qt 0
set title "Δt PDF with fitted exponential"
set xlabel "Δt"
set ylabel "PDF"
plot "waiting_time_hist_pdf.txt" using 1:2 with points lw 2 title "Simulation PDF", \
     f_pdf(x) with lines lw 2 title sprintf("λ̂ e^{-λ̂ t}, λ̂=%.4f", lambda_hat)
