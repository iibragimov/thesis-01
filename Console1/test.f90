!======<Testing find_line method>============
!Тестирование метода find_line, который написал Ренат Фаритович. 
!Изменение в адаптивном подборе шага в плоскости z.

subroutine test_find_line
use mod
external dw_dzeta, lines_test_stop_my, lines_test_stop2_my
logical lines_test_stop_my, lines_test_stop2_my
character(8) zone_name !название зон линий, для записи в файл
integer(4) i, k, zone_name_n, number_of_lines, port_z, port_zz
real(8) dir0, dl, shag, dk
complex(8) z0, zz0, zl(0:nmax), zlz(0:nmax), z, dw_dzeta, dz_dzeta, temp_zl_k

!инициализация параметров для нахождения линий тока
call init_lines_const

port_z = port
port_zz = port + 1
open(port_z, FILE='data/test_lines_z.dat')
open(port_zz, FILE='data/test_lines_zeta.dat')

!Начинаем с точки разветвления потока
z0 = cdexp(ii * g_o)
dir0 = g_o
call find_line(z0, dir0, -1, zl, k, nmax, dw_dzeta, lines_test_stop_my)
zone_name = 'first'
do i = 0, k
    zlz(i) = z(zl(i))
end do
call save_line(port_z, k, zlz, zone_name, 0)
call save_line(port_zz, k, zl, zone_name, 0)

temp_zl_k = zl(k)
number_of_lines = 20
shag = 0.02d0

!над разветлвяющейся линией
do zone_name_n = 1, number_of_lines + 5
    z0 = temp_zl_k * cdexp(-ii * shag * zone_name_n)    ! в плоскости zeta
    dir0 = -zarg(dw_dzeta(z0))
    call find_line(z0, dir0, 1, zl, k, nmax, dw_dzeta, lines_test_stop2_my)
    do i = 0, k
        zlz(i) = z(zl(i))
    end do
    zone_name = char(47 + zone_name_n)
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)    
end do
!под разветлвяющейся линией
do zone_name_n = 1, number_of_lines
    z0 = temp_zl_k * cdexp(ii * shag * zone_name_n)    ! в плоскости zeta
    dir0 = -zarg(dw_dzeta(z0))
    call find_line(z0, dir0, 1, zl, k, nmax, dw_dzeta, lines_test_stop2_my)
    do i = 0, k
        zlz(i) = z(zl(i))
    end do
    zone_name = char(47 + zone_name_n)
    call save_line(port_z, k, zlz, zone_name, zone_name_n)
    call save_line(port_zz, k, zl, zone_name, zone_name_n)    
end do

!этот кусок рисует линию тока, которая от пластинки отходит
!нужно его продумать, а еще добавить выдув струи по хорошему для второй задачи, а не для первой
z0 = c1
dir0 = -zarg(dw_dzeta(z0)) !поток вроде как по направления закрылка стекает, понять нужно почему именно -zarg(dw_dzeta(z0))
call find_line(z0, dir0, 1, zl, k, nmax, dw_dzeta, lines_test_stop2_my)
do i = 0, k
    zlz(i) = z(zl(i))
end do
zone_name = 'last_line'
call save_line(port_z, k, zlz, zone_name, 0)
call save_line(port_zz, k, zl, zone_name, 0)

call save_plastin(port_z)
call save_circle(port_zz)

close(port_z)
close(port_zz)
port = port + 2
end subroutine 