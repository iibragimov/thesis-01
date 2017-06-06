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
! метод Ньютона (находит g_b, g_p)
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

subroutine find_const
!нахождение основных параметров и констант задачи
use mod
real(8) a, b, lam
!Начальное приближение
a = pi / 3
b = pi * 3 / 4 
!Нахождение g_b, g_p
call solve_equation(b * d2, a * d2)

arg_C = (g_p - g_b) * delta / d2 - delta * pi - g_p
mod_C = l_2 / (d2*d2) * dsin(g_p * d5)**(delta - 2) * dsin(g_b * d5)**(-delta)
g_a = g_b + g_p - pi
u_inf = mod_C * v_inf
betta = alpha - arg_C
g_o = d2 * betta - 3 * pi !?
!g_o = d2 * betta + pi !-4*pi
gamma = d1 * 4 * pi * u_inf * dsin(betta)
lam = l_2 / l_1
Cy = d2 * pi * lam / (lam + d1) * dsin(pi * delta + g_p - (g_p - g_b) * d5 * delta + alpha) / dsin(g_p * d5)**(2 - delta) / dsin(g_b * d5)**delta
end subroutine

function z(zz)
!z -- отображает ед. окружность в многоугольник
!zz -- компл переменная \zeta
use mod
real(8) sgn !, x 
complex(8) z, zz, C 
C = (mod_C * cdexp(arg_C * ii)) 
!вариант с выбором ветвей
!z = C * cdexp(ii * (x + g_p + pi * sgn(x - g_p)) / d2 * (d2 - delta)) * dabs(d2 * dsin((x - g_p) / d2))**(d2 - delta) * cdexp(ii * (x + g_b + pi * sgn(x - g_b)) / d2 * delta) * dabs(d2 * dsin((x - g_b) / d2))**delta * cdexp(-ii * x)
!тупой вариант
!z = C * (zz - cdexp(ii * g_p))**(d2 - delta) * (zz - cdexp(ii * g_b))**(delta) * zz**(-d1)
!исправил, чтобы не было скачков
z = C * (zz - cdexp(ii * g_p))**(d2) * ((zz - cdexp(ii * g_b)) / (zz - cdexp(ii * g_p)))**(delta) * zz**(-d1)
end function

subroutine forma_plastinki
! записывает в файл форму пластинки
! Нужно еще разбить вычисление так чтобы все крайние и угловые точки входили g_a, g_b, g_p и т.д.
use mod
integer(4) i
complex(8) z, t 
OPEN (port, FILE = 'data/forma.dat') 
write(port, *) ' VARIABLES = "X", "Y" '
do i = 1, n
    t = z(cdexp(ii * (i - 1) * d2 * pi / (n - 1)))
    tt(i) = t
    write(port, "(F12.5,' ', F12.5)") dreal(t), dimag(t)
end do
close(port)
port = port + 1
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
OPEN (port, FILE = 'data/vs.dat') 
write(port, *) ' VARIABLES = "X", "Y" '
do i = 1, n
    !gj = (i - 1) * d2 * pi / (n - 1)
    gj = d2 * pi - (i - 1) * d2 * pi / (n - 1) 
    write(port, "(F12.5,' ', F12.5)") sg(i), v_g(gj)
end do
close(port)
port = port + 1
end subroutine

subroutine Cy_delta
! процедура построения графиков зависимости подъемной силы от дельта(угла наклона щитка) при различных длинах щитка
use mod
integer(4) i, k
k = 21
alpha = pi / 18 
OPEN (port, FILE = 'data/cy.dat') 

l_2 = 0.10d0
write(port, *) ' VARIABLES = "X", "Y" '
write(port, *) 'ZONE T = "l2 = ', l_2, '", I = ',k ,', F=POINT'
do i = 1, k
    delta = d1 + d5 * (i - 1) / (k - 1) 
    a = pi / 3
    b = pi * 3 / 4 
    call solve_equation(b * d2, a * d2)
    call find_const
    write(port, "(F12.5,' ', F12.5)") delta, Cy     
end do

l_2 = 0.20d0
write(port, *) ' VARIABLES = "X", "Y" '
write(port, *) 'ZONE T = "l2 = ', l_2, '", I = ',k ,', F=POINT'
do i = 1, k
    delta = d1 + d5 * (i - 1) / (k - 1) 
    a = pi / 3
    b = pi * 3 / 4 
    call solve_equation(b * d2, a * d2)
    call find_const
    write(port, "(F12.5,' ', F12.5)") delta, Cy     
end do

l_2 = d5
write(port, *) ' VARIABLES = "X", "Y" '
write(port, *) 'ZONE T = "l2 = ', l_2, '", I = ',k ,', F=POINT'
do i = 1, k
    delta = d1 + d5 * (i - 1) / (k - 1) 
    a = pi / 3
    b = pi * 3 / 4 
    call solve_equation(b * d2, a * d2)
    call find_const
    write(port, "(F12.5,' ', F12.5)") delta, Cy     
end do

close(port)
port = port + 1 
end subroutine

subroutine solve(x, y, z)
!главный решатель
!x -- delta угол закрылка
!y -- alpha угол атаки
!z -- l_2 длина закрылка
use mod 
real(8) x, y, z
    delta = x
    alpha = y
    l_2 = z
    call find_const
    call forma_plastinki
    call v_s
end subroutine 

!==============< Вспомогательные функции >===============!

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

!=============< Построение линий тока >==================!
!NOTICE: 
!       1. в моих обозначениях z -- это плоскость z, zz -- плоскость zeta
!       2. в обозначениях Р.Ф. наоборот
!Подключить файлы func.f90, lines_func.f90, lines_mod.f90
!Подготовить 2 функции dz_d\zeta и dw_d\zeta
!Так как функция dz_dzeta дает скачки, то попробуем мы используем функцию z(dw_dzeta(zeta))
!Подготовить 2 функции для остановки итерационного процесса
!Или написать универсальную функцию для lines_test_stop
!Записать отдельной зоной в файл саму пластинкув виде 3х точекч

function dw_dzeta(zz)
!нахождение функции dw_d\zeta
!zz -- \zeta комплексная переменная
!разделил, чтоб не таким длинным было выражение
use mod
complex(8) zz, dw_dzeta
complex(8) first, second, third
first = u_inf * cdexp(-ii * betta)
second = u_inf * cdexp(ii * betta) / zz / zz
third = gamma / (d2 * pi * ii * zz)
dw_dzeta = first - second - third
end function 

function dw_dz(zz)
!комплексно-сопряженная скорость dw_d\zeta
use mod
complex(8) dw_dzeta, dz_dzeta, dw_dz, zz
dw_dz = dw_dzeta(zz) / dz_dzeta(zz)
end function 

function dz_dzeta(zz)
!нахождение функции dz_d\zeta, поынтегральная функция
!разделил на скобки, чтоб не запупаться
!zz -- \zeta комплексная переменная
use mod
complex(8) zz, dz_dzeta
complex(8) C, zz_a, zz_p, zz_b, zz_c1
C = mod_C * cdexp(ii * arg_C)
zz_a = (zz - cdexp(ii * g_a))
zz_c1 = (zz - c1)
zz_p = (zz - cdexp(ii * g_p))
zz_b = (zz - cdexp(ii * g_b))
dz_dzeta = (C * zz_a * zz_c1 * zz_p) * (zz_b / zz_p)**(delta) / (zz_b * zz**(d2))
end function 

!zz_p = (zz - cdexp(ii * g_p))
!zz_b = (zz - cdexp(ii * g_b)) 
!zz_delta = ((zz - cdexp(ii * g_b)) / (zz - cdexp(ii * g_p)))**(delta) ! ((zeta - zeta_b) / (zeta - zeta_p))**(delta)
!dz_dzeta = C * zz_a * zz_c1 * zz_p * zz_delta / zz_b / zz**(d2)
!

!продумать нужно 
function lines_test_stop_my(z,zz)
!условие остановки при построении линии тока
!z -- что такое?
!zz -- ?
use mod
complex(8) z,zz
logical lines_test_stop_my
lines_test_stop_my = (cdabs(z) > 20.0d0) .OR. (cdabs(zz) > 20.0d0)
end

function lines_test_stop2_my(z,zz)
!условие остановки при построении линии тока
!z -- ?
!zz -- ?
use mod
complex(8) z,zz
logical lines_test_stop2_my
lines_test_stop2_my = (cdabs(z) > 20.0d0) .OR. (cdabs(zz) > 20.0d0)
end

subroutine save_line(current_port, k, array, zone_name, zone_name_n)
!запись в файл линии тока по зонам
!k -- кол-во точек линии тока
!zl -- zl(0:nmax) комплексный массив линии тока
!zone_name -- название зоны
!zone_name_n -- ?
use mod 
integer(4) i, k, zone_name_n, current_port
character(8) zone_name
complex(8) array(0:nmax)
write(current_port, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, zone_name_n, k+1
do i = 0, k
    write(current_port, "(F9.5, ' ', F9.5)") dreal(array(i)), dimag(array(i))
enddo
end subroutine 

subroutine save_plastin(port_z)
!вспомогательная процедура для сохранения пластинки в файл
use mod
integer(4) i, port_z
character(8) zone_name
complex(8) z, t 
zone_name = 'plastina'
write(port_z, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, 34, 3
write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_a))), dimag(z(cdexp(ii * g_a)))
write(port_z, "(F9.5, ' ', F9.5)") dreal(z(cdexp(ii * g_b))), dimag(z(cdexp(ii * g_b)))
write(port_z, "(F9.5, ' ', F9.5)") dreal(z(c1)), dimag(z(c1))
end subroutine 

subroutine save_circle(port_zz)
!вспомогательная процедура для сохранения окружности в файл
use mod 
integer(4) i, port_zz, k
character(8) zone_name 
real(8) gj
zone_name = 'circle'
k = 201
write(port_zz, "('ZONE T=""', A5, i2, '"", I=', i4, ', F=POINT')") zone_name, 35, 201
do i = 1, k
    gj = d2 * pi - (i - 1) * d2 * pi / (k - 1)
    write(port_zz, "(F9.5, ' ', F9.5)") dcos(gj), dsin(gj)
end do
end subroutine 

subroutine current_lines !(x, y, z)
!построение линий тока
!x -- delta угол закрылка
!y -- alpha угол атаки
!z -- l_2 длина закрылка
use mod
external dz_dzeta, dw_dzeta, lines_test_stop_my, lines_test_stop2_my
logical lines_test_stop_my, lines_test_stop2_my
character(8) zone_name !название зон линий, для записи в файл
integer(4) i, k, zone_name_n, number_of_lines, port_z, port_zz
real(8) dir0, dl, shag, dk
complex(8) z0, zz0, zl(0:nmax), zlz(0:nmax), z, dw_dzeta, dz_dzeta, temp_zl_k
!Описание переменных
!dw_dzeta -- dw_d\zeta 
!dz_dzeta -- dz_d\zeta 
!lines_test_stop, lines_test_stop2 -- функции (условие) остановки итерационного процесса, задают границы линий тока
!kol_lines -- количество линий тока
!z0 -- начальная точка построения линии тока zeta
!zz0 -- начальная точка построения линии тока в плосткости z
!dir0 -- начальное направление скорости линии тока в плосткости zeta dir0 = -zarg(dw_dzeta())
!kdir - 1 - по потоку, -1 - против потока   
!upper_bound -- верхняя граница линий тока в плоскости zeta
!lower_bound -- нижняя граница линий тока в плоскости zeta
!zl - выходной массив точек в плоскости zeta
!zlz - выходной массив точек в плоскости z
!zone_name_n -- индекс линий тока, (чтоб в файле не путать)
!z_port -- индекс для записи в файл линий тока в плоскости z
!zz_port -- индекс для записи в файл линий тока в плоскости zeta 

!нахождение параметров/констант задачи
!call find_const !возможно придется убрать и искать вне этой процедуры,тогда не нужно будет передавать xyz

!инициализация параметров для нахождения линий тока
call init_lines_const

!TODO:
!Нужно взять первую линию ту, которая врезается в пластинку и раздваивается.
!чтобы ее найти нужно найти точку раветвления потока z0 в точке О и направление для интегрирования
!Потом вызвать процедуру в противоположную сторону интегрирования
!нужно использовать find_lines2
!написать lines_test_stop универсальный

port_z = port
port_zz = port + 1
open(port_z, FILE='data/current_lines_z.dat')
open(port_zz, FILE='data/current_lines_zeta.dat')

!Предполагаемо нашли точку разветвления потока
z0 = cdexp(ii * g_o)
zz0 = z(z0)
dir0 = g_o
call find_line2(z0, zz0, dir0, -1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop_my)
zone_name = 'first'
!do i = 0, k
!    zl(i) = z(zl(i))
!end do
call save_line(port_z, k, zlz, zone_name, 0)
call save_line(port_zz, k, zl, zone_name, 0)

temp_zl_k = zl(k)
number_of_lines = 20
shag = 0.02d0
!upper_bound = dimag(zl(k)) * d2 *5   
!lower_bound = dreal(zl(k))

!над разветлвяющейся линией
do zone_name_n = 1, number_of_lines + 5
    z0 = temp_zl_k * cdexp(-ii * shag * zone_name_n)    ! в плоскости zeta
    zz0 = z(z0)                                         ! в плоскости z
    dir0 = -zarg(dw_dzeta(z0))
    call find_line2(z0, zz0, dir0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop2_my)
    zone_name = char(47 + zone_name_n)
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)    
end do
!под разветлвяющейся линией
do zone_name_n = 1, number_of_lines
    z0 = temp_zl_k * cdexp(ii * shag * zone_name_n)    ! в плоскости zeta
    zz0 = z(z0)                                         ! в плоскости z
    dir0 = -zarg(dw_dzeta(z0))
    call find_line2(z0, zz0, dir0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop2_my)
    zone_name = char(47 + zone_name_n)
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)    
end do

!этот кусок рисует линию тока, которая от пластинки отходит
!нужно его продумать, а еще добавить выдув струи по хорошему для второй задачи, а не для первой
z0 = c1
zz0 = z(c1) !вроде задняя кромка закрылка 
dir0 = -zarg(dw_dzeta(z0)) !может стоит zarg(z0) взять, поток вроде как по направления закрылка стекает, понять нужно почему именно -zarg(dw_dzeta(z0))
call find_line2(z0, zz0, dir0, 1, zl, zlz, k, nmax, dw_dzeta, dz_dzeta, lines_test_stop2_my)
zone_name = 'line' !char(47 + zone_name_n + 1)
call save_line(port_z, k, zlz, zone_name, 0)
call save_line(port_zz, k, zl, zone_name, 0)

call save_plastin(port_z)
call save_circle(port_zz)

close(port_z)
close(port_zz)
port = port + 2
end subroutine 
