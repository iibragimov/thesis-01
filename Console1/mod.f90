module mod
use gen_mod

!real(8),PARAMETER :: delta = 1.20d0           ! delta*pi - угол между щитком и пластинкой. delta = [0, 2] 
!real(8),PARAMETER ::  alpha = pi / 6       ! угол между пластинкой и горизонтом
real(8) delta, alpha, l_2


integer(4),PARAMETER :: n = 301
integer(4),PARAMETER :: nmax = 1201
real(8),PARAMETER :: l_1 = 1.0d0                ! длина пластинки
!real(8),PARAMETER :: l_2 = 0.10d0               ! длина щитка
real(8),PARAMETER :: v_inf = d1                 ! скорость на бесконечности
!real(8),PARAMETER :: eps = 0.000010d0          ! эпсюлон, беск. мал. величина

real(8) g_a, g_b, g_p, arg_C, mod_C, u_inf, gamma, betta, Cy 
real(8) vs(n), sg(n)
complex(8) tt(n) ! форма пластинки
end module mod