lambda = 1.0
n = 50

I(x) = (x > 0) ? (lambda * x - 1 - log(lambda * x)) : 1/0
Z = (exp(1)**n / lambda) * (gamma(n+1) / n**(n-1))

P(x) = (x > 0) ? exp(-n * I(x / n)) / Z : 0

plot [0:200] P(x) title "Normalized LDF PDF (n = 10)" lw 2