### theory_only.gnuplot  — underdamped free Langevin (cold start) theory outputs

reset
set decimalsign locale
set datafile separator whitespace

## ----- Physical + run params (edit as you like) -----
m     = 1.0
gamma = 0.5
kBT   = 1.0
T     = 50.0          # observation time (your T_max)

# plot / grid ranges
x_min = -50.0; x_max = 50.0
v_min =  -5.0; v_max =  5.0
NX = 101; NY = 101    # grid density for joint PDF

## ----- Derived theory moments (cold start x(0)=0, v(0)=0) -----
tau = m/gamma
a   = exp(-T/tau)
D   = kBT/gamma

sigv2 = (kBT/m) * (1.0 - a*a)
sigx2 = 2.0*D * ( T - 2.0*tau*(1.0 - a) + 0.5*tau*(1.0 - a*a) )
covxv = D * (1.0 - a)*(1.0 - a)

# covariance matrix S = [[sigx2, covxv], [covxv, sigv2]]
detS  = sigx2*sigv2 - covxv*covxv
inv11 =  sigv2/detS
inv12 = -covxv/detS
inv22 =  sigx2/detS

rho_th   = covxv / sqrt(sigx2*sigv2)
theta_th = 0.5*atan2( 2.0*covxv, sigx2 - sigv2 )   # radians
deg(x)   = x*180.0/pi

## ----- Theory PDFs -----
P_joint(x,v) = exp( -0.5*(inv11*x*x + 2.0*inv12*x*v + inv22*v*v) ) / (2.0*pi*sqrt(detS))
P_v_T(v)     = exp(-v*v/(2.0*sigv2)) / sqrt(2.0*pi*sigv2)          # finite-time velocity
P_v_MB(v)    = sqrt(m/(2.0*pi*kBT)) * exp( -m*v*v/(2.0*kBT) )      # equilibrium MB (T→∞)
P_x_T(x)     = exp(-x*x/(2.0*sigx2)) / sqrt(2.0*pi*sigx2)          # finite-time position

## ----- Write theory grid: joint_theory50.txt -----
set isosamples NX, NY
set xrange [x_min:x_max]
set yrange [v_min:v_max]
set table "joint_theory50.txt"
    splot P_joint(x,y)           # columns: x  y  P(x,y)
unset table

## ----- Write 1D marginals -----
set samples 600

# Velocity: finite-time
set table "velocity_theory50.txt"
    plot [v_min:v_max] P_v_T(x)  # columns: v  P_v_T(v)
unset table

# Velocity: Maxwell–Boltzmann (equilibrium reference)
set table "velocity_MB.txt"
    plot [v_min:v_max] P_v_MB(x) # columns: v  P_MB(v)
unset table

# Position: finite-time
set table "position_theory50.txt"
    plot [x_min:x_max] P_x_T(x)  # columns: x  P_x_T(x)
unset table

## ----- Axes through origin using theory angle -----
Lx = sqrt(sigx2)
Lv = sqrt(sigv2)
L  = (Lx>Lv)?Lx:Lv   # ~1σ; scale if you want longer/shorter axes

ct = cos(theta_th); st = sin(theta_th)
set print "axes_theta_theory50.txt"
    print "# theory axes at theta_th (x v)"
    # along theta
    print sprintf("% .12e % .12e", -L*ct, -L*st)
    print sprintf("% .12e % .12e",  L*ct,  L*st)
    print ""
    # perpendicular
    print sprintf("% .12e % .12e", -L*st,  L*ct)
    print sprintf("% .12e % .12e",  L*st, -L*ct)
unset print

## ----- Diagnostics + normalization check (coarse Riemann sum on grid) -----
# approximate integral of P_joint over the grid (should be ~1)
dx = (x_max - x_min)/(NX)
dv = (v_max - v_min)/(NY)
norm = 0.0
do for [i=0:NY-1] {
    v = v_min + (i+0.5)*dv
    do for [j=0:NX-1] {
        x = x_min + (j+0.5)*dx
        norm = norm + P_joint(x,v)*dx*dv
    }
}

set print "diagnostics_Tmax50.txt"
    print sprintf("t = %.6f", T)
    print "---- THEORY (cold start) ----"
    print sprintf("sigma_v^2(th) = %.12e", sigv2)
    print sprintf("sigma_x^2(th) = %.12e", sigx2)
    print sprintf("Cov[x,v](th)  = %.12e", covxv)
    print sprintf("rho(th)       = %.12e", rho_th)
    print sprintf("theta(th)[deg]= %.12f", deg(theta_th))
    print "---- NORMALIZATION CHECK ----"
    print sprintf("∫ P_theory dx dv ≈ %.12e", norm)
unset print

## ----- (Optional) Preview: contour-only plot of theory -----
# comment/uncomment this block to preview on screen
set term qt
set view map
unset surface
unset clabel
set contour base
set cntrparam levels discrete 0.002,0.004,0.006,0.008,0.010,0.012
splot P_joint(x,y) w l lc rgb "red" lw 1.6 t "Theory"

