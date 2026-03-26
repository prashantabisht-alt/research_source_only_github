set title 'Position Distribution (Third-Order Lambda)'
set xlabel 'Position x (units)'
set ylabel 'P(x) (units^{-1})'
set key top right
plot 'pdf_third_order_lambda.txt' using 1:2 with points title 'Simulated', \
     'pdf_third_order_lambda.txt' using 1:3 with lines title 'Theoretical'
