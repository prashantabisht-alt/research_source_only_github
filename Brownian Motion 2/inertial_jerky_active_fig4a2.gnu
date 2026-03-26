set logscale xy
set xlabel "t / τI"
set ylabel "MSD / (γ^2 v0^2 τI^6 / λ^2)"
plot "msd_active_taupratio_1e-3.txt" u 1:2 w l lw 2 t "active τp=τI/1000", \
     0.0278*x**6    w l dt 2 lw 2 t "t^6 prefactor", \
     0.0002*x**5    w l dt 2 lw 2 t "t^5 prefactor", \
     0.001333*x**3  w l dt 2 lw 2 t "t^3 prefactor"

