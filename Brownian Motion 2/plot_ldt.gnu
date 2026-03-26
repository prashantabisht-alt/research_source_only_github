# Gnuplot script to visualize Large Deviation Theory concepts



# --- Plot Aesthetics ---
set title "CLT vs. Large Deviation: Gaussian vs. Exponential Noise" font ",16"
set xlabel "Final Position (x)"
set ylabel "Probability P(x)"

# --- CRITICAL SETTING: Logarithmic y-axis to see the tails ---
set logscale y
set yrange [1e-5:1]
set grid

# --- Define a theoretical Gaussian function to fit ---
# P(x) = A * exp(-x**2 / (2*sigma**2))
# The variance of a sum of N independent variables is N * (variance of one variable).
# Var(Gaussian) = 1. Var(Sum) = Tmax * 1 = 100. So sigma = sqrt(100) = 10.
# A is a normalization constant, adjust it to match the peak of your data.
sigma = 10.0
A = 0.04
p(x) = A * exp(-x**2 / (2 * sigma**2))

# --- Plot Command ---
plot 'pdf_gaussian.txt' with points pt 7 ps 0.5 title 'Gaussian Noise', \
     'pdf_exponential.txt' with points pt 5 ps 0.5 title 'Exponential Noise', \
     p(x) with lines lw 2 lc 'black' title 'Theoretical Gaussian (CLT Prediction)'
