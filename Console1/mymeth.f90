subroutine find_const
use mod
real(8) lam
arg_C = (g_p - g_b) * delta / d2 - delta * pi - g_p
mod_C = l_2 / (d2*d2) * dsin(g_p * d5)**(delta - 2) * dsin(g_b * d5)**(-delta)
g_a = g_b + g_p - pi
u_inf = mod_C * v_inf
betta = alpha - arg_C
gamma = d1 * 4 * pi * u_inf * dsin(betta)
lam = l_2 / l_1
Cy = d2 * pi * lam / (lam + d1) * dsin(pi * delta + g_p - (g_p - g_b) * d5 * delta + alpha) / dsin(g_p * d5)**(2 - delta) / dsin(g_b * d5)**delta
end subroutine

function f1(x, y)
! первая функция, для нахождения g_p -> x, g_b -> y
use mod 
real(8) f1, x, y
f1 = dtan(x / d2) + (d2 - delta) / delta * dtan(y / d2)
end function

function f2(x, y)
! вторая функция, для нахождения g_p -> x, g_b -> y
use mod
real(8) f2, x, y 
f2 = dsin(x / d2)**(d2 - delta) * dsin(y / d2)**delta / dcos(y / d2)**(d2 - delta) / dcos(pi - x / d2)**delta - l_2 / l_1   
end function 
    
function df1_dx(x, y)
! производная от f1 по x
use mod
real(8) df1_dx, x, y
df1_dx = d5 * d1 / dcos(x / d2)**d2
end function

function df1_dy(x, y)
! производная от f1 по y
use mod
real(8) df1_dy, x, y
df1_dy = (d2 - delta) / delta / d2 / dcos(y / d2)**d2
end function

function df2_dx(x, y)
! производная от f2 по x
use mod
real(8) df2_dx, x, y
df2_dx = d5 * (-dcos(x / d2))**(-d1 - delta) * ((delta - d1) * dcos(x) - d1) * dcos(y / d2)**(delta - d2) * dsin(x / d2)**(d1 - delta) * dsin(y / d2)**delta
end function

function df2_dy(x, y)
! производная от f2 по y
use mod
real(8) df2_dy, x, y
df2_dy = d5 * (-dcos(x / d2))**(-delta) * dcos(y / d2)**(delta - 3 * d1) * (d1 + (delta - d1) * dcos(y)) * dsin(x / d2)**(d2 - delta) * dsin(y / d2)**(delta - d1)
end function

subroutine maximum(x, y, m)
! наибольшее число из х и у записывает в m, стоит переписать ее в виде ФУНКЦИИ
use mod
real(8) x, y, m
if (dabs(x) > dabs(y)) then
    m = x
else
    m = y
end if
end subroutine

subroutine solve_equation(x0, y0) !(f1, f2, df1_dx, df1_dy, df2_dx, df2_dy, x0, y0)
! метод Ньютона
use mod 
external f1, f2, df1_dx, df1_dy, df2_dx, df2_dy
real(8) f1, f2, df1_dx, df1_dy, df2_dx, df2_dy, x0, y0, a, b, c, d, opr, eps
eps = d1
do while (eps > 0.00001d0)
    a = df1_dx(x0, y0)
    b = df1_dy(x0, y0)
    c = df2_dx(x0, y0)
    d = df2_dy(x0, y0)
    opr = a * d - b * c !? нужен ли модуль
    x0 = x0 - (d * f1(x0, y0) - b * f2(x0, y0)) / opr
    y0 = y0 - (-c * f1(x0, y0) + a * f2(x0, y0)) / opr
    a = f1(x0,y0)
    b = f2(x0,y0)
    call maximum(a, b, eps)
end do

g_p = x0
g_b = y0
end subroutine

function z(x)
use mod
real(8) x, sgn 
complex(8) z, C
C = (mod_C * cdexp(arg_C * ii)) 
z = C * cdexp(ii * (x + g_p + pi * sgn(x - g_p)) / d2 * (d2 - delta)) * dabs(d2 * dsin((x - g_p) / d2))**(d2 - delta) * cdexp(ii * (x + g_b + pi * sgn(x - g_b)) / d2 * delta) * dabs(d2 * dsin((x - g_b) / d2))**delta * cdexp(-ii * x)
end function

subroutine forma_plastinki
! записывает в файл форму пластинки
use mod
integer(4) i
complex(8) z, t 
OPEN (1,FILE = 'forma.dat') 
write(1,*) ' VARIABLES = "X", "Y" '
do i = 1, n
    t = z((i - 1) * d2 * pi / (n - 1))
    tt(i) = t
    write(1,"(F12.5,' ', F12.5)") dreal(t), dimag(t)
end do
close(1)
end subroutine 

function v_g(x)
! ф. скорости V от гамма
use mod
real(8) v_g, x
v_g = v_inf * dcos(x / d2 - betta) / dabs(dcos((x - (g_p + g_b)) / d2)) * ( dabs(dsin((x - g_p) / d2) / dsin((x - g_b) / d2)))**(delta - d1)
end function

subroutine v_s
! процедура нахождения распр. скорости по пластинке и запись в файл
! нужно сделать так чтобы принимал параметры пластинки
use mod
integer(4) i
real(8) v_g, gj
sg(1) = d0
do i = 2, n   
    sg(i) = sg(i-1) + cdabs(tt(i) - tt(i-1))
end do
print *, sg(n)
OPEN (3,FILE = 'vs.dat') 
write(3,*) ' VARIABLES = "X", "Y" '
do i = 1, n
    !gj = (i - 1) * d2 * pi / (n - 1)
    gj = d2 * pi - (i - 1) * d2 * pi / (n - 1) 
    write(3,"(F12.5,' ', F12.5)") sg(i), v_g(gj)
end do
close(3)
end subroutine

subroutine Cy_delta
! процедура построения графиков зависимости подъемной силы от дельта(угла наклона щитка) при различных длинах щитка
use mod
integer(4) i, k
k = 21
alpha = pi / 18 
OPEN (2,FILE = 'cy.dat') 

l_2 = 0.10d0
write(2,*) ' VARIABLES = "X", "Y" '
write(2,*) 'ZONE T = "l2 = ', l_2, '", I = ',k ,', F=POINT'
do i = 1, k
    delta = d1 + d5 * (i - 1) / (k - 1) 
    a = pi / 3
    b = pi * 3 / 4 
    call solve_equation(b * d2, a * d2)
    call find_const
    write(2,"(F12.5,' ', F12.5)") delta, Cy     
end do

l_2 = 0.20d0
write(2,*) ' VARIABLES = "X", "Y" '
write(2,*) 'ZONE T = "l2 = ', l_2, '", I = ',k ,', F=POINT'
do i = 1, k
    delta = d1 + d5 * (i - 1) / (k - 1) 
    a = pi / 3
    b = pi * 3 / 4 
    call solve_equation(b * d2, a * d2)
    call find_const
    write(2,"(F12.5,' ', F12.5)") delta, Cy     
end do

l_2 = d5
write(2,*) ' VARIABLES = "X", "Y" '
write(2,*) 'ZONE T = "l2 = ', l_2, '", I = ',k ,', F=POINT'
do i = 1, k
    delta = d1 + d5 * (i - 1) / (k - 1) 
    a = pi / 3
    b = pi * 3 / 4 
    call solve_equation(b * d2, a * d2)
    call find_const
    write(2,"(F12.5,' ', F12.5)") delta, Cy     
end do

close(2)
end subroutine


!===========================< * >===================================!

function sred(x1,x2)
!среднее арифметическое
use mod
real(8) sred, x1, x2
sred = (x1 + x2) / d2
end function

function sgu(x)
!сгущает точки к середине
use mod
real(8) x, sgu
sgu = (x - d5) * (x - d5) * (x - d5) / (d5 * d5) + d5
!sgu=d1 - ( atan( tan(( (1-x)*(1-x) )*pi/d2) ) )*d2/pi
end function

function sgn(x)
use mod
real(8) sgn, x
if (x < d0) then
    sgn = -d1
else
    sgn = d1
end if 
end function

subroutine solve(x, y, z)
use mod 
real(8) a, b, x, y, z
    delta = x
    alpha = y
    l_2 = z
    a = pi / 3
    b = pi * 3 / 4 
    call solve_equation(b * d2, a * d2)
    call find_const
    call forma_plastinki
    call v_s
end subroutine 
    
