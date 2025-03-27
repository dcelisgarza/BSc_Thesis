! For numerical eigenvalue calculations http://people.sc.fsu.edu/~jburkardt/f_src/eispack/eispack.html

program main
implicit none
!call window_test
!call kdelta_test
!call print_PES
call sofar
end program main

subroutine sofar
integer, parameter :: n=4
integer f, i
double precision w(n), avgnum(n/2), denom
call print_pes(n)
call MC_Average(i, f, w, avgnum, denom, n)
end subroutine sofar

subroutine print_PES(n)
! Testing testing all PES
implicit none
integer n, flag
double precision a, b, c, d, eo, r, h(n), e(n/2), omega, dh(n), de(n)
do flag = 1, n-1
   if (flag == 1) then
      r = -4.0
      a = 0.01
      b = 1.6
      c = 0.005
      d = 1
      open(unit=1, file='sac.dat')
      open(unit=2, file='dsac.dat')
      call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
      call deriv_PES(flag, a, b, c, d, eo, r, dh, de, n)
      write(1,*) r, h(1), h(2), h(3), h(4), e(1), e(2)
      write(2,*) r, dh(1), dh(2), dh(3), dh(4)
      do while(r<4)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
         call deriv_PES(flag, a, b, c, d, eo, r, dh, de, n)
         write(1,*) r, h(1), h(2), h(3), h(4), e(1), e(2)
         write(2,*) r, dh(1), dh(2), dh(3), dh(4)
      end do
   else if (flag == 2) then
      r  = -8.0
      a  = 0.1
      b  = 0.28
      c  = 0.015
      d  = 0.06
      eo = 0.05
      open(unit=1, file='dac.dat')
      open(unit=2, file='ddac.dat')
      call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
      call deriv_PES(flag, a, b, c, d, eo, r, dh, de, n)
      write(1,*) r, h(1), h(2), h(3), h(4), e(1), e(2)
      write(2,*) r, dh(1), dh(2), dh(3), dh(4)
      do while(r<8)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
         call deriv_PES(flag, a, b, c, d, eo, r, dh, de, n)
         write(1,*) r, h(1), h(2), h(3), h(4), e(1), e(2)
         write(2,*) r, dh(1), dh(2), dh(3), dh(4)
      end do
   else if (flag == 3) then
      r = -10.0
      a = 6e-4
      b = 0.1
      c = 0.9
      open(unit=1, file='ec.dat')
      open(unit=2, file='dec.dat')
      call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
      call deriv_PES(flag, a, b, c, d, eo, r, dh, de, n)
      write(1,*) r, h(1)*50, h(2), h(3), h(4)*50, e(1), e(2)
      write(2,*) r, dh(1), dh(2), de(1), de(2)
      do while(r<10)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
         call deriv_PES(flag, a, b, c, d, eo, r, dh, de, n)
         write(1,*) r, h(1)*50, h(2), h(3), h(4)*50, e(1), e(2)
         write(2,*) r, dh(1), dh(2), de(1), de(2)
      end do
   end if
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt') 
102 format(5(1pe12.3))
end subroutine print_PES

subroutine sc_motion(t,x,dx,n)
!==========================================================!
! Equations of motion of single avoided crossing           !
! Daniel Celis Garza 21-22 Jan 2015                        !
!----------------------------------------------------------!
! x(1) = R     dx(1) = dR/dt                               !
! x(2) = P     dx(2) = dP/dt                               !
! x(3) = x1    dx(3) = dx1/dt                              !
! x(4) = p1    dx(4) = dp1/dt                              !
! x(5) = x2    dx(5) = dx2/dt                              !
! x(6) = p2    dx(6) = dp2/dt                              !
!==========================================================!
implicit none
integer n
double precision t, x(n+2), dx(n+2), h(n), dh(n), e(n), de(n), a, b, c, d, eo, omega, mu
a    = 0.01
b    = 1.6
c    = 0.005
d    = 1.0
call PES(1, a, b, c, d, eo, x(1), h, e, omega, n)
call deriv_PES(1, a, b, c, d, eo, x(1), dh, de, n)
! Speed adjustment routine goes here.
dx(1) = x(2)/mu
dx(2) = -((dh(1)+dh(4))/2 + (x(3)+x(4)-x(5)-x(6))*(dh(1)-dh(4))/4 + dh(3)*(x(3)*x(5) + x(4)*x(6)))
dx(3) = x(4)*(h(1)-h(4))/2 + x(6)*h(3)
dx(4) = -(x(3)*(h(1)-h(4))/2 + x(5)*h(3))
dx(5) = -x(6)*(h(1)-h(4))/2 + x(4)*h(3)
dx(6) = -(-x(5)*(h(1)-h(4))/2 + x(3)*h(3))
end subroutine sc_motion

subroutine dc_motion(t,x,dx,n)
!==========================================================!
! Equations of motion of single avoided crossing           !
! Daniel Celis Garza 21-22 Jan 2015                        !
!----------------------------------------------------------!
! x(1) = R     dx(1) = dR/dt                               !
! x(2) = P     dx(2) = dP/dt                               !
! x(3) = x1    dx(3) = dx1/dt                              !
! x(4) = p1    dx(4) = dp1/dt                              !
! x(5) = x2    dx(5) = dx2/dt                              !
! x(6) = p2    dx(6) = dp2/dt                              !
!==========================================================!
implicit none
integer n
double precision t, x(n+2), dx(n+2), h(n), dh(n), e(n), de(n), a, b, c, d, eo, omega, mu
a  = 0.1
b  = 0.28
c  = 0.015
d  = 0.06
eo = 0.05
call PES(2, a, b, c, d, eo, x(1), h, e, omega, n)
call deriv_PES(2, a, b, c, d, eo, x(1), dh, de, n)
! Speed adjustment routine goes here.
dx(1) = x(2)/mu
dx(2) = -((dh(1)+dh(4))/2 + (x(3)+x(4)-x(5)-x(6))*(dh(1)-dh(4))/4 + dh(3)*(x(3)*x(5) + x(4)*x(6)))
dx(3) = x(4)*(h(1)-h(4))/2 + x(6)*h(3)
dx(4) = -(x(3)*(h(1)-h(4))/2 + x(5)*h(3))
dx(5) = -x(6)*(h(1)-h(4))/2 + x(4)*h(3)
dx(6) = -(-x(5)*(h(1)-h(4))/2 + x(3)*h(3))
end subroutine dc_motion

subroutine ec_motion(t,x,dx,n)
!==========================================================!
! Equations of motion of single avoided crossing           !
! Daniel Celis Garza 21-22 Jan 2015                        !
!----------------------------------------------------------!
! x(1) = R     dx(1) = dR/dt                               !
! x(2) = P     dx(2) = dP/dt                               !
! x(3) = x1    dx(3) = dx1/dt                              !
! x(4) = p1    dx(4) = dp1/dt                              !
! x(5) = x2    dx(5) = dx2/dt                              !
! x(6) = p2    dx(6) = dp2/dt                              !
!==========================================================!
implicit none
integer n
double precision t, x(n+2), dx(n+2), h(n), dh(n), e(n), de(n), a, b, c, d, eo, omega, mu, c1, c2, dc1
a  = 6e-4
b  = 0.1
c  = 0.1
mu = 2000
call PES(3, a, b, c, d, eo, x(1), h, e, omega, n)
call deriv_PES(3, a, b, c, d, eo, x(1), dh, de, n)
! Speed adjustment routine goes here (i think).
c1    = (x(2)+omega*(x(4)*x(5)-x(6)*x(3)))/mu
c2    = e(1)-e(2)
dc1   = (dh(2)*(h(1)-h(4))+h(2)*(dh(4)-dh(1)))*(x(4)*x(5)-x(6)*x(3))/(4*h(2)**2+(h(1)-h(4))**2)
dx(1) = c1
dx(2) = -c1*dc1-(de(1)-de(2))*(x(3)**2+x(4)**2-x(5)**2-x(6)**2)/4-(de(1)+de(2))/2
dx(3) = c1*omega*x(5)+x(4)*c2/2
dx(4) = c1*omega*x(6)-x(3)*c2/2
dx(5) = -c1*omega*x(3)-x(6)*c2/2
dx(6) = -c1*omega*x(4)+x(5)*c2/2
end subroutine ec_motion

subroutine deriv_PES(flag, a, b, c, d, eo, r, dh, de, n)
!===============================================================!
! Analytical derivatives of the PES used in Hamilton's eqs      !
!---------------------------------------------------------------!
! Daniel Celis Garza 21 Jan 2015                                !
!===============================================================!
implicit none
integer n, flag
double precision a, b, c, d, eo, r, dh(n), de(n/2)
if (flag == 1) then
   if (r > 0) then
      dh(1) = a*b*exp(-b*r) 
   else
      dh(1) = a*b*exp(b*r)
   end if
   dh(2) = -2*c*d*r*exp(-d*r**2)
   dh(3) = dh(2)
   dh(4) = -dh(1)
else if (flag == 2) then
   dh(1) = 0.0
   dh(2) = -2*c*d*r*exp(-d*r**2)
   dh(3) = dh(2)
   dh(4) = 2*a*b*r*exp(-b*r**2)
else if (flag == 3) then
   dh(1) = 0
   dh(4) = 0
   if (r > 0) then
      de(1) = -b**2*c*exp(-2*c*r)*(2*exp(c*r)-1)/sqrt(a**2+b**2*(2-exp(-c*r))**2)
      dh(2) = b*c*exp(-c*r)
   else
      de(1) = -b**2*c*exp(2*c*r)/sqrt(a**2+b**2*exp(2*c*r))
      dh(2) = b*c*exp(c*r)
   end if
   de(2) = -de(1)
   dh(3) = dh(2)
end if
end subroutine deriv_PES

subroutine PES(flag, a, b, c, d, eo, r, h, e, omega, n)
!=========================================================!
! Potential energy surfaces as defined by Miller          !
!---------------------------------------------------------!
! Daniel Celis Garza 16 Jan 2015, 20 Jan 2015             !
!---------------------------------------------------------!
! flag = 1  ->  single avoided crossing (diabatic)        !
! flag = 2  ->  dual avoided crossing (diabatic)          !
! flag = 3  ->  extended coupling (adiabatic)             !
! flag = 4  ->  spin-boson model condensed phase dynamics !
!---------------------------------------------------------!
! h(n)      -> diabatic matrix elements                   !
! e(n/2)    -> adiabatic matrix elements (eigenvalues of  !
!              the diabatic matrix)                       !
!---------------------------------------------------------!
! h(1) = H11    h(2) = H12    h(3) = H21    h(4) = H22    !
! e(1) = E1     e(2) = E2                                 !
!=========================================================!
implicit none
integer n, flag
double precision a, b, c, d, eo, r, h(n), e(n/2), omega
! Single Avoided Crossing
if (flag == 1) then
   if (r > 0) then
      h(1) = a*(1-exp(-b*r))
      e(1) = -exp(-d*r**2)*sqrt(c**2+exp(2*d*r**2)*(a*(1-exp(-b*r)))**2)
      e(2) = -e(1)
   else
      h(1) = -a*(1-exp(b*r))
      e(1) = -exp(-d*r**2)*sqrt(c**2+exp(2*d*r**2)*(-a*(1-exp(b*r)))**2)
      e(2) = -e(1)
   end if
   h(2) = c*exp(-d*r**2)
   h(3) = h(2)
   h(4) = -h(1)
! Double Avoided Crossing
else if (flag == 2) then
   h(1) = 0
   h(2) = c*exp(-d*r**2)
   h(3) = h(2)
   h(4) = -a*exp(-b*r**2) + eo
   e(1) = exp(-(b+d)*r**2)*(-a*exp(d*r**2)+exp((b+d)*r**2)*eo-sqrt(4*c**2*exp(2*b*r**2)+exp(2*d*r**2)*(a-exp(b*r**2)*eo)**2))/2
   e(2) = exp(-(b+d)*r**2)*(-a*exp(d*r**2)+exp((b+d)*r**2)*eo+sqrt(4*c**2*exp(2*b*r**2)+exp(2*d*r**2)*(a-exp(b*r**2)*eo)**2))/2
! Extended Coupling
else if (flag == 3) then
   h(1) = -a
   if (r > 0) then
      h(2) = b*(2-exp(-c*r))
   else
      h(2) = b*exp(c*r)
   end if
   h(3) = h(2)
   h(4) = -h(1)
   omega = atan(2*h(2)/(h(1)-h(4)))/2
   if (r > 0) then
      e(1) = -sqrt(a**2 + (b*(2-exp(-c*r)))**2)
      e(2) = -e(1)
   else
      e(1) = -sqrt(a**2 + (b*exp(c*r))**2)
      e(2) = -e(1)
   end if
end if
end subroutine PES

subroutine MC_Average(i, f, w, avgnum, denom, n)
!=======================================================!
! Window function as defined on equation 11 of          !
! meyer & miller's  paper (base paper for the thesis)   !
! W_N (n) = 1/delta_n * h(delta_n/2 - |n-N|)            !
!-------------------------------------------------------!
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
! Monte-Carlo Averaging procedure
! Daniel Celis Garza 20 Jan 2015
! w(1) to w(n/2) we have the window functions for the initial conditions
! w(n/2 + 1) to w(n) we have the window functions for the final conditions
implicit none
integer kdelta, n, i, f, k, run, eqn
double precision heaviside, arg, del_ea, ea, gamma, rndnum, w(n), num(n/2), denom, avgnum(n/2)
num(1) = 1.0
num(2) = 1.0
avgnum(1) = 0.0
avgnum(2) = 0.0
denom  = 0.0
i = 1
gamma  = (sqrt(3.0)-1.0)/2.0
del_ea = 2.0*gamma
call init_random_seed()
do f = 1, n/2
   do run = 1, 100 ! Monte-Carlo averaging.
      do k = 1, n/2
         ! Calculating all the window functions for the initial electronic state.
         call random_number(rndnum)
         eqn  = kdelta(i, k)
         ea   = eqn + gamma*(2.0*rndnum-1.0) ! initial ea in paper it's eq 13
         arg  = del_ea/2.0 - abs(ea-eqn)
         w(k) =  heaviside(arg)/del_ea
         ! Call the subroutine which calculates the trajectory.
         ! Calculating all the window functions for the final electronic state.
         call random_number(rndnum)
         eqn  = kdelta(f, n/2+k)
         ea   = eqn + gamma*(2.0*rndnum-1.0) ! this should be the final ea, eq (14) in the paper
         arg  = del_ea/2.0 - abs(ea-eqn)
         w(n/2+k) =  heaviside(arg)/del_ea
         ! Calculating the numerator before monte carlo averaging.
         num(f) = num(f)*w(k)*w(n/2+k)
      end do
      ! Adding the numerator for monte-carlo averaging.
      avgnum(f) = avgnum(f)+num(f)
      num(f) = 1.0
   end do
   avgnum(f) = avgnum(f)/100
   denom     = denom + avgnum(f) ! Obtaining the denominator with the MC-averaged numerator
end do
end subroutine MC_Average

!subroutine trans_prob
!=======================================================!
! Monte Carlo Averaging as defined on equation 15 of    !
! meyer & miller's paper (base paper for the thesis)    !
!-------------------------------------------------------!
! Daniel Celis Garza 15 Jan 2015 - 20 Jan 2015          !
!-------------------------------------------------------!
!implicit none
!integer arg(2), f, i, k
!double precision kdelta
!i      = 1                     ! Initial electronic state.
!f      = 2                     ! Final electronic state.
!gamma  = (sqrt(3.0)-1.0)/2.0   ! gamma parameter.
!del_ea = 2.0*gamma             ! Delta_n
!eqn    = kdelta(arg(1), arg(2))
!ea     = eqn + gamma*(2.0*rndnum-1.0)
!h_arg  = del_ea/2.0 - abs(ea-eqn)
!end subroutine trans_prob

function heaviside(arg)
!=======================================================!
! Heaviside function                                    !
! Daniel Celis Garza 15 Jan 2015                        !
!=======================================================!
implicit none
double precision heaviside, arg
heaviside = 0.0
if (arg >= 0.0) heaviside = 1.0
end function heaviside

function kdelta(i,j)
!=======================================================!
! Kroencker delta function                              !
! Daniel Celis Garza 15 Jan 2015                        !
! delta_ij                                              !
!=======================================================!
implicit none
integer i, j, kdelta
kdelta = 0
if (i == j) kdelta = 1
end function kdelta

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

!subroutine kdelta_test
! 15 Jan 2015
!implicit none
!integer kdelta, arg(2)
!arg(1) = 2.5*2
!arg(2) = 5
!print *, kdelta(arg(1), arg(2))
!end subroutine kdelta_test
