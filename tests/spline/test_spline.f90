program test_spline

implicit none

integer :: n !Number of data points
real(8), dimension(:), allocatable :: t, y, ypp !Data abscissa, ordinate, 2nd deriv
integer :: nval !Number of interpolation points
real(8), dimension(:), allocatable :: tval, yval, ypval, yppval
real(8) :: dt
integer :: fu
integer :: i

!Read in data
open(newunit=fu, file='input.txt', action='read', status='old')

!Ignore the first line
read(fu, *)

!Get number of data points
read(fu, *) n

allocate( t(n), y(n), ypp(n) )
t = 0.d0; y = 0.d0; ypp = 0.d0

do i = 1, n
    read(fu, *) t(i), y(i)
end do

close(fu)

nval = 124
allocate( tval(nval), yval(nval), ypval(nval), yppval(nval) )
tval = 0.d0; yval = 0.d0; ypval = 0.d0; yppval = 0.d0
tval(1) = minval(t); tval(nval) = maxval(t)
dt = ( tval(nval) - tval(1) )/(nval-1)
do i = 2, (nval-1)
    tval(i) = tval(1) + (i-1)*dt
end do

!Construct interpolating cubic spline with natural BCs
call spline_cubic_set(n, t, y, 2, 0.d0, 2, 0.d0, ypp)

!Evaluate spline
do i = 1, nval
    call spline_cubic_val(n, t, y, ypp, tval(i), yval(i), ypval(i), yppval(i))
end do

!Write output
open(newunit=fu, file='csp.txt', action='write', status='replace')
write(fu, *) '#t  y  yp  ypp'
do i = 1, nval
    write(fu, '(*(g0.12,2x))') tval(i), yval(i), ypval(i), yppval(i)
end do
close(fu)

!Linear spline
yval = 0.d0; ypval = 0.d0
do i = 1, nval
    call spline_linear_val(n, t, y, tval(i), yval(i), ypval(i))
end do

!Write output for linear spline
open(newunit=fu, file='lsp.txt', action='write', status='replace')
write(fu, *) '#t  y  yp'
do i = 1, nval
    write(fu, '(*(g0.12,2x))') tval(i), yval(i), ypval(i)
end do
close(fu)

end program test_spline
