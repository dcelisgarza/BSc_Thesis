! For numerical eigenvalue calculations http://people.sc.fsu.edu/~jburkardt/f_src/eispack/eispack.html

module parameters
double precision :: w, kd
end module

program main
implicit none
!call window_test
!call kdelta_test
call test_PES
end program main

subroutine test_PES
! Testing testing all PES
implicit none
integer, parameter :: n = 4
integer flag
double precision a, b, c, d, eo, r, h(n)
flag = 2
a    = 6e-4
b    = 0.1
c    = 0.09
d    = 0.06
eo   = 0.05
r    = -10.0
call PES(flag, a, b, c, d, eo, r, h, n)

open(unit=1, file = 'extended_coup.dat')
! Write initial values.
write(1,*) r, h(1)*50.0, h(2), h(3), h(4)*50.0
do while(r < 8)
   ! Write the values for other r's.
   r = r + 0.1
   call PES(flag, a, b, c, d, eo, r, h, n)
   write(1,*) r, h(1)*50.0, h(2), h(3), h(4)*50.0
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt') 
102 format(5(1pe12.3))
end subroutine test_PES

subroutine PES(flag, a, b, c, d, eo, r, h, n)
!=========================================================!
! Potential energy surfaces as defined by Miller          !
!---------------------------------------------------------!
! Daniel Celis Garza 16 Jan 2015                          !
!---------------------------------------------------------!
! flag = 0  ->  single avoided crossing                   !
! flag = 1  ->  dual avoided crossing                     !
! flag = 2  ->  extended coupling                         !
! flag = 3  ->  spin-boson model condensed phase dynamics !
! h(1) = H11    h(2) = H12    h(3) = H21    h(4) = H22    !
!=========================================================!
implicit none
integer n, flag
double precision a, b, c, d, eo, r, h(n)
if (flag == 0) then
   if (r > 0) then
      h(1) = a*(1-exp(-b*r))
   else
      h(1) = -a*(1-exp(b*r))
   end if
   h(2) = c*exp(-d*r*r)
   h(3) = h(2)
   h(4) = -h(1)
else if (flag == 1) then
   h(1) = 0
   h(2) = c*exp(-d*r*r)
   h(3) = h(2)
   h(4) = -a*exp(-b*r*r) + eo
else if (flag == 2) then
   h(1) = -a
   if (r > 0) then
      h(2) = b*(2-exp(-c*r))
   else
      h(2) = b*exp(c*r)
   end if
   h(3) = h(2)
   h(4) = -h(1)
end if
end subroutine PES

!subroutine trans_prob
!==========================================================================================!
! Calculate transition probabilities from initial to final electronic state                !
! P_{f <- i} = \frac{\prod_{k=1}^{F} W_{\delta_{fk}} \dot \prod_{k=1}^{F} W_{\delta_{ik}}} !
!    {\sum_{k=1}^{F} \prod_{k=1}^{F} W_{\delta_{fk}} \dot \prod_{k=1}^{F} W_{\delta_{ik}}} !
!------------------------------------------------------------------------------------------!
! Daniel Celis Garza 15 Jan 2015 - 
!==========================================================================================! 
!implicit none
!use parameters
!double precision f, i, k
!end subroutine trans_prob

subroutine window_test
! Test subroutine for the window
! 15 Jan 2015
use parameters
implicit none
double precision h_arg, del_ea, ea, eqn, gamma, rndnum
call init_random_seed()
call random_number(rndnum)
gamma  = (sqrt(3.0)-1.0)/2.0
del_ea = 2.0*gamma
eqn    = 1.0
ea     = eqn + gamma*(2.0*rndnum-1.0)
h_arg  = del_ea/2.0 - abs(ea-eqn)
call window_func(h_arg,del_ea)
print *, h_arg, del_ea/2, -abs(ea-eqn), w
print *, ea, eqn - del_ea/2.0, eqn + del_ea/2.0
end subroutine window_test

subroutine window_func(arg,del_ea)
use parameters
!=======================================================!
! Window function as defined on equation 11 of          !
! meyer & miller's  paper (base paper for the thesis)   !
! W_N (n) = 1/delta_n * h(delta_n/2 - |n-N|)            !

! Daniel Celis Garza 15 Jan 2015                        !
!-------------------------------------------------------!
! delta_n -> del_ea (delta electronic action variable)  !
! n       -> ea     (electronic action variable)        !
! N       -> eqn    (electronic quantum number)         !
!-------------------------------------------------------!
! For there to be a histogram bin the following         !
! condition must be met:                                !
! eqn - delta_n/2 <= ea <= eqn - del_ea                 !
!=======================================================!
implicit none
double precision heaviside, arg, del_ea, ea, eqn
w = heaviside(arg)/del_ea
end subroutine window_func

function heaviside(arg)
!=======================================================!
! Heaviside function                                    !
! Daniel Celis Garza 15 Jan 2015                        !
!=======================================================!
implicit none
double precision heaviside, arg
heaviside = 0.0
if (arg .ge. 0.0) heaviside = 1.0
end function heaviside

subroutine kdelta_test
! 15 Jan 2015
use parameters
implicit none
double precision kdelta, arg(2)
arg(1) = -0.999999975
arg(2) = -0.999999975
kd = kdelta(arg(1), arg(2))
print *, kd
end subroutine kdelta_test

function kdelta(arg1, arg2)
!=======================================================!
! Kroencker delta function                              !
! Daniel Celis Garza 15 Jan 2015                        !
! arg(1) = f                                            !
! arg(2) = i                                            !
!=======================================================!
implicit none
double precision kdelta, arg1, arg2
kdelta = 0.0
if (arg1 .eq. arg2) kdelta = 1.0
end function kdelta

subroutine motion(t,x,dx,n)
use parameters
implicit none
integer n
double precision t, x(n), dx(n)
end subroutine motion

subroutine rk4n(f,ti,tf,xi,xf,n)
!=======================================================!
! Solve n first-order ODEs or n/2 second-order.         !
! Runge-Kutta 4                                         !
! Daniel Celis Garza 10 Dec. 2014                       !
!-------------------------------------------------------!
! f(t,x,dx,n) = parametrised equations of motion        !
! ti          = independent variable @ step i           !
! tf          = independent variable @ step i+1         !
! xi()        = array of dependent variables @ step i   !
! xf()        = array of dependent variables @ step i+1 !
! n           = # of parametrised equations of motion   !
!-------------------------------------------------------!
! j           = counter variable                        !
! h           = time step                               !
! t           = current time                            !
! ki()        = array of Runge-Kutta k's                !
!=======================================================!
implicit none
integer n, j
double precision ti, tf, h, t, xi(n), xf(n), k1(n), k2(n), k3(n), k4(n), x(n), dx(n)

h = tf - ti ! This will be replaced by the error-correcting routine.
t = ti

! Calculate k1
call f(t,xi,dx,n)            ! Call the equations of motion
do j = 1, n                  ! Go through all equations of motion
   k1(j) = h*dx(j)           ! Calculate the value of k1 for each equation of motion
   x(j)  = xi(j) + k1(j)/2.0 ! Calculate the next value of x for each equation (to be used in the next ki)
end do

! Calculate k2
call f(t+h/2.0,x,dx,n)
do j = 1, n
   k2(j) = h*dx(j)
   x(j)  = xi(j) + k2(j)/2.0
end do

! Calculate k3
call f(t+h/2.0,x,dx,n)
do j = 1, n
   k3(j) = h*dx(j)
   x(j)  = xi(j) + k3(j)/2.0
end do

! Calculate k4 and xf
call f(t+h,x,dx,n)
do j = 1, n
   k4(j) = h*dx(j)
   xf(j)  = xi(j) + k1(j)/6.0 + k2(j)/3.0 + k3(j)/3.0 + k4(j)/6.0
end do
end subroutine rk4n

subroutine init_random_seed()
!=======================================================!
!-----------Standard F95 Random seed routine------------!
! rndnum can be any n x m array or a scalar             !
! REAL :: rndnum(n,m)                                   !
! CALL init_random_seed()                               !
! CALL RANDOM_NUMBER(rndnum)                            !
!=======================================================!
use iso_fortran_env, only: int64
implicit none
integer, allocatable :: seed(:)
integer :: i, n, un, istat, dt(8), pid
integer(int64) :: t

call random_seed(size = n)
allocate(seed(n))
! First try if the OS provides a random number generator
open(newunit=un, file="/dev/urandom", access="stream", &
 form="unformatted", action="read", status="old", iostat=istat)
if (istat == 0) then
   read(un) seed
   close(un)
else
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
   call system_clock(t)
   if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
   end if
   pid = getpid()
   t = ieor(t, int(pid, kind(t)))
   do i = 1, n
      seed(i) = lcg(t)
   end do
end if
call random_seed(put=seed)
contains
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
function lcg(s)
integer :: lcg
integer(int64) :: s
if (s == 0) then
s = 104729
else
s = mod(s, 4294967296_int64)
end if
s = mod(s * 279470273_int64, 4294967291_int64)
lcg = int(mod(s, int(huge(0), int64)), kind(0))
end function lcg
end subroutine init_random_seed
