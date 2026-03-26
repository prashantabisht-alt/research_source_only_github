program ldf_calculator
  implicit none
  integer, parameter :: nbins = 2000
  real*8, dimension(nbins) :: x, p_gauss, p_exp
  real*8 :: xmin, xmax, bin_width, t, D, b_theory
  integer :: i, ios, n_tail_g, n_tail_e
  real*8 :: slope_g, intercept_g, slope_e, intercept_e
  real*8, dimension(:), allocatable :: x2_tail_g, logp_tail_g, x_tail_e, logp_tail_e
  real*8 :: sum_x, sum_logp, sum_x2, sum_x_logp
  real*8 :: I_theory_g, I_theory_e, s

  ! Parameters
  t = 10.0d0    ! Total time (Tmax * dt = 10000 * 0.001)
  D = 1.0d0     ! Diffusion coefficient
  xmin = -100.0d0
  xmax = 100.0d0
  bin_width = (xmax - xmin) / dble(nbins)
  b_theory = sqrt(D / 0.001d0)  ! Theoretical b from variance match, b = sqrt(D / dt)

  ! Initialize arrays
  x = 0.0d0
  p_gauss = 0.0d0
  p_exp = 0.0d0

  ! Read data files
  open(unit=10, file='pdf_gaussian.txt', status='old', iostat=ios)
  if (ios /= 0) stop 'Error opening pdf_gaussian.txt'
  do i = 1, nbins
    read(10, *, iostat=ios) x(i), p_gauss(i)
    if (ios /= 0) exit
  end do
  close(10)

  open(unit=11, file='pdf_exponential.txt', status='old', iostat=ios)
  if (ios /= 0) stop 'Error opening pdf_exponential.txt'
  do i = 1, nbins
    read(11, *, iostat=ios) x(i), p_exp(i)
    if (ios /= 0) exit
  end do
  close(11)

  ! Count tail points (|x| > 10)
  n_tail_g = 0
  n_tail_e = 0
  do i = 1, nbins
    if (abs(x(i)) > 10.0d0 .and. p_gauss(i) > 1.0d-20) n_tail_g = n_tail_g + 1
    if (abs(x(i)) > 10.0d0 .and. p_exp(i) > 1.0d-20) n_tail_e = n_tail_e + 1
  end do
  allocate(x2_tail_g(n_tail_g), logp_tail_g(n_tail_g))
  allocate(x_tail_e(n_tail_e), logp_tail_e(n_tail_e))

  ! Extract tail data for Gaussian
  n_tail_g = 0
  do i = 1, nbins
    if (abs(x(i)) > 10.0d0 .and. p_gauss(i) > 1.0d-20) then
      n_tail_g = n_tail_g + 1
      x2_tail_g(n_tail_g) = x(i)**2
      logp_tail_g(n_tail_g) = log(p_gauss(i))
    end if
  end do

  ! Linear fit for Gaussian (ln P vs x^2)
  sum_x = 0.0d0
  sum_logp = 0.0d0
  sum_x2 = 0.0d0
  sum_x_logp = 0.0d0
  do i = 1, n_tail_g
    sum_x = sum_x + x2_tail_g(i)
    sum_logp = sum_logp + logp_tail_g(i)
    sum_x2 = sum_x2 + x2_tail_g(i)**2
    sum_x_logp = sum_x_logp + x2_tail_g(i) * logp_tail_g(i)
  end do
  slope_g = (n_tail_g * sum_x_logp - sum_x * sum_logp) / (n_tail_g * sum_x2 - sum_x**2)
  intercept_g = (sum_logp - slope_g * sum_x) / dble(n_tail_g)

  ! Extract tail data for Exponential
  n_tail_e = 0
  do i = 1, nbins
    if (abs(x(i)) > 10.0d0 .and. p_exp(i) > 1.0d-20) then
      n_tail_e = n_tail_e + 1
      x_tail_e(n_tail_e) = abs(x(i))
      logp_tail_e(n_tail_e) = log(p_exp(i))
    end if
  end do

  ! Linear fit for Exponential (ln P vs |x|)
  sum_x = 0.0d0
  sum_logp = 0.0d0
  sum_x2 = 0.0d0
  sum_x_logp = 0.0d0
  do i = 1, n_tail_e
    sum_x = sum_x + x_tail_e(i)
    sum_logp = sum_logp + logp_tail_e(i)
    sum_x2 = sum_x2 + x_tail_e(i)**2
    sum_x_logp = sum_x_logp + x_tail_e(i) * logp_tail_e(i)
  end do
  slope_e = (n_tail_e * sum_x_logp - sum_x * sum_logp) / (n_tail_e * sum_x2 - sum_x**2)
  intercept_e = (sum_logp - slope_e * sum_x) / dble(n_tail_e)

  ! Compare with theoretical LDF at a sample s
  s = 1.0d0  ! Example s = x/t, e.g., x=10 at t=10
  I_theory_g = (s**2) / (4.0d0 * D)  ! Theoretical Gaussian LDF
  I_theory_e = abs(s) / b_theory     ! Theoretical Exponential LDF

  ! Output results with comparison
  print *, 'Large Deviation Function Analysis (t = ', t, '):'
  print *, 'Gaussian Noise:'
  print *, '  Fit: ln P ≈ ', intercept_g, ' - ', -slope_g, ' * x^2'
  print *, '  Computed I(s) ≈ ', -slope_g * t, ' * s^2'
  print *, '  Theoretical I(s) = s^2 / (4D) = ', I_theory_g, ' at s = ', s
  print *, '  Difference: ', (-slope_g * t) - I_theory_g
  print *, 'Exponential Noise:'
  print *, '  Fit: ln P ≈ ', intercept_e, ' - ', -slope_e, ' * |x|'
  print *, '  Computed I(s) ≈ ', -slope_e * t, ' * |s|'
  print *, '  Theoretical I(s) = |s| / b ≈ ', I_theory_e, ' at s = ', s
  print *, '  Difference: ', (-slope_e * t) - I_theory_e

  deallocate(x2_tail_g, logp_tail_g, x_tail_e, logp_tail_e)

end program ldf_calculator