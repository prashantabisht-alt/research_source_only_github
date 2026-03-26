### ------------------------------------------------------------
### plot_ou.gp — Single PDF summary of OU verification outputs
### ---------------------------------------------
set grid
set key outside top center horizontal
set tics nomirror

### Output format: PDF (landscape)
set term pdfcairo size 11.0,8.5 enhanced font ",11"

### Layout: 2 rows × 2 columns
set multiplot layout 2,2 title "Ornstein–Uhlenbeck Process Summary"

# ------------------------------------------------------------
# Panel 1: PDF vs Gaussian (linear scale)
# ------------------------------------------------------------
set key left top
set title "PDF: Histogram vs Gaussian"
set xlabel "u"
set ylabel "Probability density"
unset logscale y
plot \
    "ou_hist_vs_gauss.txt" using 1:2 with lines lw 2 title "P_sim(u)", \
    "" using 1:3 with lines lw 2 dt 2 title "Gaussian theory"

# ------------------------------------------------------------
# Panel 2: Eq. (6) covariance (semilog-y)
# ------------------------------------------------------------
set key right top
set title "Eq. (6) covariance"
set xlabel "τ"
set ylabel "Covariance"
set logscale y
plot \
    "ou_cov_eq6.txt" using 2:3 with lines lw 2 title "C_emp(τ)", \
    "" using 2:4 with lines lw 2 dt 2 title "v0^2 e^{-τ/τp}"
unset logscale y

# ------------------------------------------------------------
# Panel 3: Normalized ACF (linear scale)
# ------------------------------------------------------------
set key right top
set title "Normalized ACF"
set xlabel "τ"
set ylabel "ACF"
unset logscale y
plot \
    "ou_acf_eq6_norm.txt" using 2:3 with lines lw 2 title "ACF_emp_norm", \
    "" using 2:4 with lines lw 2 dt 2 title "exp(-τ/τp)"

# ------------------------------------------------------------
# Panel 4: Time series window
# ------------------------------------------------------------
set key left top
set title "u(t) — First 5k samples"
set xlabel "t"
set ylabel "u(t)"
plot "ou_series.txt" using 2:3 every ::1::5000 with lines lw 1.5 title "u(t)"

unset multiplot
unset output
set term qt   # back to interactive mode if supported
