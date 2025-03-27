program main
  use functions, only : init_random_seed
  implicit none
  call ec_data_collection
end program main

subroutine ec_data_collection
integer, parameter :: n = 1, es = 2 ! n = nucleus, es = electronic states
! avgnumt, avgnumr = average numerator for transmission and reflection respectively
! denom = denominator; ip = initial momentum
double precision avgnumt(es), avgnumr(es), denom, ip
! Time stamps
double precision days, hours, minutes, seconds, start, finish
CALL CPU_TIME(start)
!open(unit=1, file='ec_prob.dat')
!write(1,*) "# h = 1/(0.10125*ip), mc_steps = 15000, gamma = 0.366, initial state = 2"
!write(1,*) "# Initial Momentum  ", "R 1<-2  ", "R 2<-2  ", "T 1<-2  ", "T 2<-2  "

ip = 13.0

!Initial ones.
call ec_mc_avg(avgnumr,avgnumt,denom,ip,n,es)
!write(1,*) ip, avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom

!do while(ip < 30.0)
!   ip = ip + 1.0
!   call ec_mc_avg(avgnumr,avgnumt,denom,ip,n,es)
!   write(1,*) ip, avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom
!end do

!CALL CPU_TIME(finish)
!seconds = finish-start
!days = floor(seconds/86400)
!hours = floor(seconds/3600.0)
!minutes = floor(mod(seconds/60.0,60.0))
!seconds = mod(seconds,60.0)
!write(1,*) '# Execution Time:'
!write(1,*) '# days = ', days, 'hrs = ', hours, 'min = ', minutes, 's = ', seconds
end subroutine ec_data_collection

subroutine ec_mc_avg(avgnumr,avgnumt,denom,ip,n,es)
implicit none
integer n, es ! Parameters
integer i, run, maxrun ! Counters
double precision numt(es), numr(es), avgnumt(es), avgnumr(es), denom, ip
! Initialising variables
denom   = 0.0
do concurrent (i = 1: es)
   avgnumt(i) = 0.0
   avgnumr(i) = 0.0
end do

! Parameters
maxrun = 1
do concurrent (run = 1: maxrun)
   call ec(ip, numt, numr, n, es)
   ! Summing numerators. Only those corresponding to the final state will change.
   do concurrent (i = 1: es)
      avgnumt(i) = avgnumt(i) + numt(i)
      avgnumr(i) = avgnumr(i) + numr(i)
   end do
end do

do concurrent (i = 1: es)
   avgnumt(i) = avgnumt(i)/maxrun
   avgnumr(i) = avgnumr(i)/maxrun
   denom      = denom + avgnumt(i) + avgnumr(i)
end do
end subroutine ec_mc_avg

subroutine ec(ip, numt, numr, n, es)
use functions, only : heaviside, kdelta, rk4gn!, init_random_seed
implicit none
double precision, parameter :: pi = 4.0*atan(1.0)
integer n, es ! Parameters
integer i, j ! Counters.
integer reflection ! Flags 1 = yes, 0 = no.
integer statei, statef, k ! E state related.
double precision ti, tf, dt, dumrmin, rmin, rmax, gamma, dsn ! Adjustment parameters.
double precision bni(es), bnf(es), sni(es), snf(es), q(es), rnd(es) ! E state related
double precision numt(es), numr(es), arg, wi(es), wf(es) ! MC Average.
! There are two coordinates per nucleus.
! There are two coordinates per electronic state.
! xi() = system's generalised coordinates.
double precision xi(2*(n+es)), xf(2*(n+es)), ip ! Variables.
double precision numdum
external ec_eq ! Single crossing Hamilton's eqs.

100 continue ! We return after an unsuccessful attempt.
! Initialising variables.
reflection = 0
gamma      = 0.366!(sqrt(3.0)-1.0)/2.0
dsn        = 2.0*gamma
arg        = 0.0
j          = 0
statef     = 0
numdum     = 1.0
do concurrent (k = 1: es)
   numt(k) = 0.0
   numr(k) = 0.0
   wi(k)   = 0.0
   wf(k)   = 0.0
end do

! Parameters.
ti      = 0.0
dt      = 1.0/(5.0*ip)
rmin    = -9.0
rmax    = 6.0
statei  = 2
dumrmin = rmin

! Action variables.
!call init_random_seed()
do k = 1, es
   call random_number(rnd)
   bni(k) = kdelta(statei,k)
   sni(k) = bni(k) + gamma*(2.0*rnd(1)-1.0)
   q(k)   = 2.0*pi*rnd(2)
   arg    = gamma-abs(sni(k)-bni(k))
   wi(k)  = heaviside(arg)/dsn
end do

! Initial cartesian variables
xi(1) = rmin ! Position
xi(2) = ip   ! Momentum
xi(3) = sqrt(2.0*(sni(1)+gamma))*cos(q(1))  ! x1
xi(4) = -sqrt(2.0*(sni(1)+gamma))*sin(q(1)) ! p1
xi(5) = sqrt(2.0*(sni(2)+gamma))*cos(q(2))  ! x2
xi(6) = -sqrt(2.0*(sni(2)+gamma))*sin(q(2)) ! p2

! Integrate hamilton's equations
open(unit=1,file='ec_traj_.dat')
write(1,*) ti, xi(1), xi(2), xi(3)**2.0/2.0 + xi(4)**2.0/2.0, xi(5)**2.0/2.0 + xi(6)**2.0/2.0
do while(dumrmin<=rmax)
   tf = ti + dt
   call rk4gn(ec_eq,ti,tf,xi,xf,2*(n+es))
   ! Check if the particle was reflected back.
   dumrmin = xf(1)
!print*, xi(1)-xf(1)
   !if(xf(1)-xi(1) .lt. 0.0) then
   if(dumrmin .lt. rmin) then
      reflection = 1
      exit
   end if
   ti = tf
   do concurrent (i = 1: 2*(n+es))
      xi(i) = xf(i)
   end do
   write(1,*) ti, xi(1), xi(2), xi(3)**2.0/2.0 + xi(4)**2.0/2.0, xi(5)**2.0/2.0 + xi(6)**2.0/2.0
end do

! Calculating nk(t).
do concurrent (k = 1: es)
   snf(k) = 0.5*xf(es+k+j)**2.0 + 0.5*xf(es+k+1+j)**2.0
   j      = j + 1
end do

if(1.0 .le. snf(1) .and. snf(1) .le. 1.0 + dsn .and. 0.0 .le. snf(2) .and. snf(2) .le. dsn) then
   statef  = 1
else if(0.0 .le. snf(1) .and. snf(1) .le. dsn .and. 1.0 .le. snf(2) .and. snf(2) .le. 1.0 + dsn) then
   statef  = 2
else
   ! Well, go to 100 (the start of the subroutine) seemed to do the trick...
   CLOSE ( 1, STATUS='DELETE', IOSTAT=I )
   go to 100
end if

do k = 1, es
   bnf(k) = kdelta(statef,k)
   arg    = gamma - abs(snf(k)-bnf(k))
   wf(k)  = heaviside(arg)/dsn
end do

if(reflection .eq. 0) then
   do concurrent (k = 1: es)
      numdum = numdum*wi(k)*wf(k)
   end do
   ! Only the numerator which corresponds to the final state will change, the other will remain equal to 0.0
   numt(statef) = numdum
   numdum       = 1.0
! Check for reflection.
else
   do concurrent (k = 1: es)
      numdum = numdum*wi(k)*wf(k)
   end do
   ! Only the numerator which corresponds to the final state will change, the other will remain equal to 0.0
   numr(statef) = numdum
   numdum       = 1.0
end if
end subroutine ec

subroutine ec_eq(t,x,dx,n)
!==========================================================!
! Equations of motion of single avoided crossing           !
! Daniel Celis Garza 2015                                  !
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
double precision t, x(n), dx(n), a, b, c, mu
double precision c0, c1, c2, c3, c4, c5, c6, mu4, pp!, dc1
a  = 6e-4
b  = 0.1
c  = 0.9
mu = 2000.0
mu4 = mu*4.0
pp = 2.0*x(2)

if (x(1) .ge. 0.0) then
   c0 = b*exp(-c*x(1))
   c1 = b*2.0-c0
   c2 = atan(c1/a)
   c3 = a**2.0+c1**2.0
   c4 = sqrt(c3)
   c5 = x(6)*x(3)-x(4)*x(5)
   c6 = pp+c5*c2
   dx(1) = (x(2)+0.5*c5*c2)/mu
   dx(2) = 0.5*c*c0*(-a*c5*c6/(mu*2.0*c3)+(x(3)**2.0+x(4)**2.0-x(5)**2.0-x(6)**2.0)*c1/c4)
   dx(3) = -x(5)*c2*c6/mu4-x(4)*c4
   dx(4) = -x(6)*c2*c6/mu4+x(3)*c4
   dx(5) = x(3)*c2*c6/mu4+x(6)*c4
   dx(6) = x(4)*c2*c6/mu4-x(5)*c4
else if (x(1) .lt. 0.0) then
   c1 = b*exp(c*x(1))
   c2 = atan(c1/a)
   c3 = a**2.0 + c1**2.0
   c4 = sqrt(c3)
   c5 = x(6)*x(3)-x(4)*x(5)
   c6 = pp+c5*c2
   dx(1) = (x(2)+0.5*c5*c2)/mu
   dx(2) = 0.5*c*c1*(-a*c5*c6/(mu*2.0*c3)+(x(3)**2.0+x(4)**2.0-x(5)**2.0-x(6)**2.0)*c1/c4)
   dx(3) = -x(5)*c2*c6/mu4-x(4)*c4
   dx(4) = -x(6)*c2*c6/mu4+x(3)*c4
   dx(5) = x(3)*c2*c6/mu4+x(6)*c4
   dx(6) = x(4)*c2*c6/mu4-x(5)*c4
end if
end subroutine ec_eq
