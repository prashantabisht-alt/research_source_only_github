program precision_demo
  implicit none
  ! ask the compiler: give me a real kind with at least
  ! 15 decimal digits of precision and exponent range ~10^±300
  integer, parameter :: dp = selected_real_kind(15,300)

  real(dp) :: x, y, z

  ! assign a literal with the _dp suffix so it's stored at double precision
  x = 1.0_dp
  y = 1.0_dp + 1.0e-16_dp   ! this difference is only visible in double precision
  z = y - x

  print *, "dp kind number = ", dp
  print *, "epsilon(1.0_dp) = ", epsilon(1.0_dp)
  print *, "tiny(1.0_dp)    = ", tiny(1.0_dp)
  print *, "huge(1.0_dp)    = ", huge(1.0_dp)

  print *, "x = ", x
  print *, "y = ", y
  print *, "z = y - x = ", z
end program precision_demo