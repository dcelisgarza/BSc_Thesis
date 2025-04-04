program main
 !use mpi_f08
    implicit none
    !integer :: rank, size, len
    !character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: version

    !call MPI_INIT()
    !call MPI_COMM_RANK(MPI_COMM_WORLD, rank)
    !call MPI_COMM_SIZE(MPI_COMM_WORLD, size)
    !call MPI_GET_LIBRARY_VERSION(version, len)

    !call sc_data_collection
    call dc_data_collection
    !call print_pes(n)
    !call MPI_FINALIZE()
end program main

subroutine sc_data_collection
integer, parameter :: n = 4
double precision xi(n+2), avgnumt(n/2), avgnumr(n/2), denom, ip, start, finish
double precision days, hours, minutes, seconds
CALL CPU_TIME(start)
open(unit=1, file='sc_prob.dat')
write(1,*) "# h = 1/ip, mc_steps = 15000, gamma = 0.366"
write(1,*) "# Initial Momentum  ", "R 1<-1  ", "R 2<-1  ", "T 1<-1  ", "T 2<-1  "

ip = 1.0

!Initial ones.
call sc_mc_avg(avgnumr,avgnumt,denom,ip,n)
write(1,*) ip, avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom

do while(ip < 30.0)
   ip = ip + 1.0
   call sc_mc_avg(avgnumr,avgnumt,denom,ip,n)
   write(1,*) ip, avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom
end do

CALL CPU_TIME(finish)
seconds = finish-start
days = floor(seconds/86400)
hours = floor(seconds/3600.0)
minutes = floor(mod(seconds/60.0,60.0))
seconds = mod(seconds,60.0)
write(1,*) '# Execution Time:'
write(1,*) '# days = ', days, 'hrs = ', hrs, 'min = ', minutes, 's = ', seconds
end subroutine sc_data_collection

subroutine sc_mc_avg(avgnumr,avgnumt,denom,ip,n)
implicit none
integer n, i, run, maxrun
double precision numt(n/2), numr(n/2), avgnumt(n/2), avgnumr(n/2), denom, ip
double precision xi(n+2)
! Initialising variables
denom   = 0.0
do concurrent (i = 1: n/2)
   avgnumt(i) = 0.0
   avgnumr(i) = 0.0
end do

! Parameters
maxrun = 15000
do concurrent (run = 1: maxrun)
   call sc(ip, numt, numr, n)
   ! Summing numerators. Only those corresponding to the final state will change.
   do concurrent (i = 1: n/2)
      avgnumt(i) = avgnumt(i) + numt(i)
      avgnumr(i) = avgnumr(i) + numr(i)
   end do
end do

do concurrent (i = 1: n/2)
   avgnumt(i) = avgnumt(i)/maxrun
   avgnumr(i) = avgnumr(i)/maxrun
   denom      = denom + avgnumt(i) + avgnumr(i)
end do
end subroutine sc_mc_avg

subroutine sc(ip, numt, numr, n)
implicit none
double precision, parameter :: pi = 4.0*atan(1.0)
integer n ! F x F matrix.
integer i, j ! Counters.
integer reflection, outofrange ! Flags 1 = yes, 0 = no.
integer statei, statef, k ! E state related.
double precision ti, tf, dt, dumrmin, rmin, rmax, gamma, dsn ! Adjustment parameters.
double precision kdelta, bni(n/2), bnf(n/2), sni(n/2), snf(n/2), q(n/2), rnd(n/2) ! E state related
double precision numt(n/2), numr(n/2), heaviside, arg, wi(n/2), wf(n/2) ! MC Average.
double precision xi(n+2), xf(n+2), ip ! Variables.
double precision numdum
external sc_eq ! Single crossing Hamilton's eqs.

100 continue ! We return after an unsuccessful attempt.
! Initialising variables.
reflection = 0
outofrange = 0
gamma      = 0.366!(sqrt(3.0)-1.0)/2.0
dsn        = 2.0*gamma
arg        = 0.0
j          = 0
statef     = 0
numdum     = 1.0
do concurrent (i = 1: n/2)
   numt(i) = 0.0
   numr(i) = 0.0
   wi(i)   = 0.0
   wf(i)   = 0.0
end do

! Parameters.
ti      = 0.0
dt      = 1.0/(ip)
rmin    = -4.0
rmax    = 4.0
statei  = 1
dumrmin = rmin

! Action variables.
call init_random_seed()
do k = 1, n/2
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
do while(dumrmin<=rmax)
   tf = ti + dt
   call rk4gn(sc_eq,ti,tf,xi,xf,n+2)
   ! Check if the particle was reflected back.
   dumrmin = xf(1)
   if(dumrmin .le. rmin) then
      reflection = 1
      exit
   end if
   ti = tf
   do concurrent (i = 1: n+2)
      xi(i) = xf(i)
   end do
end do

! Calculating nk(t).
do concurrent (k = 1: n/2)
   snf(k) = 0.5*xf(n/2+k+j)**2.0 + 0.5*xf(n/2+k+1+j)**2.0
   j      = j + 1
end do

if(1.0 .le. snf(1) .and. snf(1) .le. 1.0 + dsn .and. 0.0 .le. snf(2) .and. snf(2) .le. dsn) then
   statef  = 1
else if(0.0 .le. snf(1) .and. snf(1) .le. dsn .and. 1.0 .le. snf(2) .and. snf(2) .le. 1.0 + dsn) then
   statef  = 2
else
   ! Well, go to 100 (the start of the subroutine) seemed to do the trick...
   go to 100
end if

do k = 1, n/2
   bnf(k) = kdelta(statef,k)
   arg    = gamma - abs(snf(k)-bnf(k))
   wf(k)  = heaviside(arg)/dsn
end do

if(reflection .eq. 0) then
   do concurrent (k = 1: n/2)
      numdum = numdum*wi(k)*wf(k)
   end do
   ! Only the numerator which corresponds to the final state will change, the other will remain equal to 0.0
   numt(statef) = numdum
   numdum       = 1.0
! Check for reflection.
else
   do concurrent (k = 1: n/2)
      numdum = numdum*wi(k)*wf(k)
   end do
   ! Only the numerator which corresponds to the final state will change, the other will remain equal to 0.0
   numr(statef) = numdum
   numdum       = 1.0
end if
end subroutine sc

subroutine sc_eq(t,x,dx,n)
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
double precision c1, c2
a  = 0.01
b  = 1.6
c  = 0.005
d  = 1.0
mu = 2000.0
c1 = c*exp(-d*x(1)**2.0) 

dx(1) = x(2)/mu
if(x(1) .ge. 0.0) then
   c2 = a - a*exp(-b*x(1))
   dx(2) = 2.0*c1*d*x(1)*(x(4)*x(6)+x(3)*x(5))+0.5*(-x(4)**2.0+x(6)**2.0-x(3)**2.0+x(5)**2.0)*a*b*exp(-b*x(1))
   dx(3) = c1*x(6) + x(4)*c2
   dx(4) = -c1*x(5) - x(3)*c2
   dx(5) = c1*x(4) - x(6)*c2
   dx(6) = -c1*x(3) + x(5)*c2
else if(x(1) .lt. 0.0) then
   c2 = a*exp(b*x(1))-a
   dx(2) = 2.0*c1*d*x(1)*(x(4)*x(6)+x(3)*x(5))+0.5*(-x(4)**2.0+x(6)**2.0-x(3)**2.0+x(5)**2.0)*a*b*exp(b*x(1))
   dx(3) = c1*x(6) + x(4)*c2
   dx(4) = -c1*x(5) - x(3)*c2
   dx(5) = c1*x(4) - x(6)*c2
   dx(6) = -c1*x(3) + x(5)*c2
end if
!call PES(1, a, b, c, d, eo, x(1), h, e, omega, n)
!call dPES(1, a, b, c, d, eo, x(1), dh, de, n)
!dx(1) = x(2)/mu
!dx(2) = -(x(4)*x(6)+x(3)*x(6))*dh(2)+1/2*(-dh(1)-dh(4))-1/2*(x(4)**2/2-x(6)**2/2+x(3)**2/2-x(6)**2/2)*(dh(1)-dh(4))
!dx(3) = x(6)*h(2)+1/2*x(4)*(h(1)-h(4))
!dx(4) = -x(6)*h(2)- 1/2*x(3)*(h(1)-h(4))
!dx(5) = x(4)*h(2)- 1/2*x(6)*(h(1)-h(4))
!dx(6) = -x(3)*h(2)+ 1/2*x(6)*(h(1)-h(4))
end subroutine sc_eq

subroutine dc_data_collection
integer, parameter :: n = 4
double precision xi(n+2), avgnumt(n/2), avgnumr(n/2), denom, ip, start, finish, avgenergy
double precision days, hours, minutes, seconds
double precision adjuster
CALL CPU_TIME(start)
open(unit=1, file='dc_prob.dat')
write(1,*) "# h = 1/5ip, mc_steps = 15000, gamma = 0.366"
write(1,*) "# Average Initial Energy ", "R 1<-1  ", "R 2<-1  ", "T 1<-1  ", "T 2<-1  "

ip = 8.0
adjuster = 1.0
!Initial ones.
call dc_mc_avg(avgnumr,avgnumt,denom,ip,avgenergy,n)
write(1,*) log(avgenergy), avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom
!write(1,*) log(ip**2.0/4000.0), avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom

do while(ip < 105.0)
   ip = ip + 2.0*adjuster
   call dc_mc_avg(avgnumr,avgnumt,denom,ip,avgenergy,n)
   write(1,*) log(avgenergy), avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom
   !write(1,*) log(ip**2.0/4000.0), avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom
   adjuster = adjuster + 0.08
end do

CALL CPU_TIME(finish)
seconds = finish-start
days = floor(seconds/86400)
hours = floor(seconds/3600.0)
minutes = floor(mod(seconds/60.0,60.0))
seconds = mod(seconds,60.0)
write(1,*) '# Execution Time:'
write(1,*) '# days = ', days, 'hrs = ', hrs, 'min = ', minutes, 's = ', seconds
end subroutine dc_data_collection

subroutine dc_mc_avg(avgnumr,avgnumt,denom,ip,avgenergy,n)
implicit none
integer n, i, run, maxrun
double precision numt(n/2), numr(n/2), avgnumt(n/2), avgnumr(n/2), denom, ip, energy, avgenergy
double precision xi(n+2)
! Initialising variables
denom     = 0.0
avgenergy = 0.0
energy    = 0.0
do concurrent (i = 1: n/2)
   avgnumt(i) = 0.0
   avgnumr(i) = 0.0
end do

! Parameters
maxrun = 100
do concurrent (run = 1: maxrun)
   call dc(ip, numt, numr, energy, n)
   ! Summing numerators. Only those corresponding to the final state will change.
   do concurrent (i = 1: n/2)
      avgnumt(i) = avgnumt(i) + numt(i)
      avgnumr(i) = avgnumr(i) + numr(i)
   end do
   avgenergy = avgenergy + energy
end do
avgenergy = avgenergy/maxrun

do concurrent (i = 1: n/2)
   avgnumt(i) = avgnumt(i)/maxrun
   avgnumr(i) = avgnumr(i)/maxrun
   denom      = denom + avgnumt(i) + avgnumr(i)
end do
end subroutine dc_mc_avg

subroutine dc(ip, numt, numr, energy, n)
implicit none
double precision, parameter :: pi = 4.0*atan(1.0)
integer n ! F x F matrix.
integer i, j ! Counters.
integer reflection, outofrange, e_counter ! Flags 1 = yes, 0 = no.
integer statei, statef, k ! E state related.
double precision ti, tf, dt, dumrmin, rmin, rmax, gamma, dsn ! Adjustment parameters.
double precision kdelta, bni(n/2), bnf(n/2), sni(n/2), snf(n/2), q(n/2), rnd(n/2) ! E state related
double precision numt(n/2), numr(n/2), heaviside, arg, wi(n/2), wf(n/2) ! MC Average.
double precision xi(n+2), xf(n+2), ip, energy ! Variables.
double precision numdum
external dc_eq ! Single crossing Hamilton's eqs.
e_counter = 0
100 continue ! We return after an unsuccessful attempt.
! Initialising variables.
reflection = 0
outofrange = 0
gamma      = 0.366!(sqrt(3.0)-1.0)/2.0
dsn        = 2.0*gamma
arg        = 0.0
j          = 0
statef     = 0
numdum     = 1.0
do concurrent (i = 1: n/2)
   numt(i) = 0.0
   numr(i) = 0.0
   wi(i)   = 0.0
   wf(i)   = 0.0
end do

! Parameters.
ti      = 0.0
dt      = 1.0/(0.1*ip)
rmin    = -8.0
rmax    = 8.0
statei  = 1
dumrmin = rmin

! Action variables.
call init_random_seed()
do k = 1, n/2
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

energy = xi(2)**2.0/4000.0 + 0.5*(-0.1*exp(-0.28*xi(1)**2.0)+0.05)&
        +0.25*(xi(3)**2.0+xi(4)**2.0-xi(5)**2.0-xi(6)**2.0)&
        *(0.1*exp(-0.28*xi(1)**2.0)-0.05)+(xi(4)*xi(6)+xi(3)*xi(5))&
        *0.015*exp(-0.06*xi(1)**2.0)

! Integrate hamilton's equations
do while(dumrmin<=rmax)
   tf = ti + dt
   call rk4gn(dc_eq,ti,tf,xi,xf,n+2)
   ! Check if the particle was reflected back.
   dumrmin = xf(1)
   if(dumrmin .le. rmin) then
      reflection = 1
      exit
   end if
   ti = tf
   do concurrent (i = 1: n+2)
      xi(i) = xf(i)
   end do
end do

! Calculating nk(t).
do concurrent (k = 1: n/2)
   snf(k) = 0.5*xf(n/2+k+j)**2.0 + 0.5*xf(n/2+k+1+j)**2.0
   j      = j + 1
end do

if(1.0 .le. snf(1) .and. snf(1) .le. 1.0 + dsn .and. 0.0 .le. snf(2) .and. snf(2) .le. dsn) then
   statef  = 1
else if(0.0 .le. snf(1) .and. snf(1) .le. dsn .and. 1.0 .le. snf(2) .and. snf(2) .le. 1.0 + dsn) then
   statef  = 2
else
   ! Well, go to 100 (the start of the subroutine) seemed to do the trick...
   go to 100
end if

do k = 1, n/2
   bnf(k) = kdelta(statef,k)
   arg    = gamma - abs(snf(k)-bnf(k))
   wf(k)  = heaviside(arg)/dsn
end do

if(reflection .eq. 0) then
   do concurrent (k = 1: n/2)
      numdum = numdum*wi(k)*wf(k)
   end do
   ! Only the numerator which corresponds to the final state will change, the other will remain equal to 0.0
   numt(statef) = numdum
   numdum       = 1.0
! Check for reflection.
else
   do concurrent (k = 1: n/2)
      numdum = numdum*wi(k)*wf(k)
   end do
   ! Only the numerator which corresponds to the final state will change, the other will remain equal to 0.0
   numr(statef) = numdum
   numdum       = 1.0
end if

!energy = xf(2)**2.0/4000.0 + 0.5*(-0.1*exp(-0.28*xf(1)**2.0)+0.05)&
!        +0.25*(xf(3)**2.0+xf(4)**2.0-xf(5)**2.0-xf(6)**2.0)&
!        *(0.1*exp(-0.28*xf(1)**2.0)-0.05)+(xf(4)*xf(6)+xf(3)*xf(5))&
!        *0.015*exp(-0.06*xf(1)**2.0)
end subroutine dc

subroutine dc_eq(t,x,dx,n)
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
double precision c1, c2
a  = 0.1
b  = 0.28
c  = 0.015
d  = 0.06
eo = 0.05
mu = 2000.0
c1 = 0.5*(a*exp(-b*x(1)**2.0)-eo)
c2 = c*exp(-d*x(1)**2.0)

dx(1) = x(2)/mu
dx(2) = 0.5*exp(-(b+d)*x(1)**2.0)*x(1)*(4*c*d*exp(b*x(1)**2.0)*(x(4)*x(6)+x(3)*x(5))&
        +a*b*exp(d*x(1)**2.0)*(-2.0+x(4)**2.0-x(6)**2.0+x(3)**2.0-x(5)**2.0))
dx(3) = c1*x(4)+c2*x(6)
dx(4) = -c1*x(3)-c2*x(5)
dx(5) = c2*x(4)-c1*x(6)
dx(6) = -c2*x(3)+c1*x(5)
!call PES(2, a, b, c, d, eo, x(1), h, e, omega, n)
!call dPES(2, a, b, c, d, eo, x(1), dh, de, n)
!dx(1) = x(2)/mu
!dx(2) = -((dh(1)+dh(4))/2 + (x(3)+x(4)-x(5)-x(6))*(dh(1)-dh(4))/4 + dh(3)*(x(3)*x(5) + x(4)*x(6)))
!dx(3) = x(4)*(h(1)-h(4))/2 + x(6)*h(3)
!dx(4) = -(x(3)*(h(1)-h(4))/2 + x(5)*h(3))
!dx(5) = -x(6)*(h(1)-h(4))/2 + x(4)*h(3)
!dx(6) = -(-x(5)*(h(1)-h(4))/2 + x(3)*h(3))
end subroutine dc_eq

subroutine ec_eq(t,x,dx,n)
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
double precision c0, c1, c2, c3, c4
a  = 6e-4
b  = 0.1
c  = 0.1
mu = 2000.0

if (x(1) .ge. 0.0) then
   c0 = b*exp(-c*x(1))
   c1 = b*2.0 - c0
   c2 = atan(c1/a)
   c3 = a**2.0+c1**2.0
   c4 = sqrt(c3)
   dx(1) = (x(2)+0.5*(x(6)*x(3)-x(4)*x(5))*c2)/mu
   dx(2) = 0.5*c0*c*((a*(x(6)*x(3)-x(4)*x(5))*(-2.0*x(2)+(x(4)*x(5)-x(6)*x(3))*c2))&
        /(2.0*mu*c3)+((x(4)**2.0-x(6)**2.0+x(3)**2.0-x(5)**2.0)*c1)/c4)
   dx(3) = 1/(4.0*mu)*x(5)*c2*(-2.0*x(2)+(x(4)*x(5)-x(6)*x(3))*c2)-x(4)*c4
   dx(4) = -(1/(4.0*mu))*x(6)*c2*(2.0*x(2)+(x(6)*x(3)-x(4)*x(5))*c2)+x(3)*c4
   dx(5) = 1/(4.0*mu)*x(3)*c2*(2.0*x(2)+(x(6)*x(3)-x(4)*x(5))*c2)+x(6)*c4
   dx(6) = 1/(4.0*mu)*x(4)*c2*(2.0*x(2)+(x(6)*x(3)-x(4)*x(5))*c2)-x(5)*c4
else if(x(1) .lt. 0.0) then
   c1 = b*exp(c*x(1))
   c2 = atan(c1/a)
   c3 = a**2.0+c1**2.0
   c4 = sqrt(c3)
   dx(1) = (x(2)+0.5*(x(6)*x(3)-x(4)*x(5))*c2)/mu
   dx(2) = 0.5*c1*c*((a*(x(6)*x(3)-x(4)*x(5))*(-2.0*x(2)+(-x(6)*x(3)+ x(4)*x(5))*c2))&
        /(2.0*mu*c3)+((x(4)**2.0-x(6)**2.0+x(3)**2.0-x(5)**2.0)*c1)/c4)
   dx(3) = 1/(4.0*mu)*x(5)*c2*(-2.0*x(2)+(x(4)*x(5)-x(6)*x(3))*c2)-x(4)*c4
   dx(4) = -(1/(4.0*mu))*x(6)*c2*(2.0*x(2)+(x(6)*x(3)-x(4)*x(5))*c2)+x(3)*c4
   dx(5) = 1/(4.0*mu)*x(3)*c2*(2.0*x(2)+(x(6)*x(3)-x(4)*x(5))*c2)+x(6)*c4
   dx(6) = 1/(4.0*mu)*x(4)*c2*(2.0*x(2)+(x(6)*x(3)-x(4)*x(5))*c2)-x(5)*c4
end if


!call PES(3, a, b, c, d, eo, x(1), h, e, omega, n)
!call dPES(3, a, b, c, d, eo, x(1), dh, de, n)
!c1    = (x(2)+omega*(x(4)*x(5)-x(6)*x(3)))/mu
!c2    = e(1)-e(2)
!dc1   = (dh(2)*(h(1)-h(4))+h(2)*(dh(4)-dh(1)))*(x(4)*x(5)-x(6)*x(3))/(4*h(2)**2+(h(1)-h(4))**2)
!dx(1) = c1
!dx(2) = -c1*dc1-(de(1)-de(2))*(x(3)**2+x(4)**2-x(5)**2-x(6)**2)/4-(de(1)+de(2))/2
!dx(3) = c1*omega*x(5)+x(4)*c2/2
!dx(4) = c1*omega*x(6)-x(3)*c2/2
!dx(5) = -c1*omega*x(3)-x(6)*c2/2
!dx(6) = -c1*omega*x(4)+x(5)*c2/2
end subroutine ec_eq

function heaviside(x)
!=======================================================!
! Heaviside function                                    !
! Daniel Celis Garza 15 Jan 2015                        !
!=======================================================!
implicit none
double precision heaviside, x
heaviside = 0.0
if (x .ge. 0.0) heaviside = 1.0
end function heaviside

function kdelta(i,j)
!=======================================================!
! Kroencker delta function                              !
! Daniel Celis Garza 15 Jan 2015                        !
! delta_ij                                              !
!=======================================================!
implicit none
integer i, j
double precision kdelta
kdelta = 0.0
if (i .eq. j) kdelta = 1.0
end function kdelta

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

subroutine dPES(flag, a, b, c, d, eo, r, dh, de, n)
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
end subroutine dPES

subroutine rk4gn(f,ti,tf,xi,xf,n)
!=======================================================!
! Solve n first-order ODEs or n/2 second-order.         !
! Runge-Kutta 4 Gill's method                           !
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
double precision ti, tf, h, t, xi(n), xf(n), k1(n), k2(n), k3(n), k4(n), x(n), dx(n), sqr2

h = tf - ti ! This will be replaced by the error-correcting routine.
t = ti
sqr2 = sqrt(2.0)
! Calculate k1
call f(t,xi,dx,n)            ! Call the equations of motion
do concurrent (j = 1: n)                  ! Go through all equations of motion
   k1(j) = h*dx(j)           ! Calculate the value of k1 for each equation of motion
   x(j)  = xi(j) + k1(j)/2.0 ! Calculate the next value of x for each equation (to be used in the next ki)
end do

! Calculate k2
call f(t+h/2.0,x,dx,n)
do concurrent (j = 1: n)
   k2(j) = h*dx(j)
   x(j)  = xi(j) + (sqr2-1.0)*k1(j)/2.0 + (1.0-sqr2/2.0)*k2(j)
end do

! Calculate k3
call f(t+h/2.0,x,dx,n)
do concurrent (j = 1: n)
   k3(j) = h*dx(j)
   x(j)  = xi(j) - sqr2*k2(j)/2.0 + (1.0+sqr2/2.0)*k3(j)
end do

! Calculate k4 and xf
call f(t+h,x,dx,n)
do concurrent (j = 1: n)
   k4(j) = h*dx(j)
   xf(j)  = xi(j) + (k1(j) + (2.0-sqr2)*k2(j) + (2+sqr2)*k3(j) + k4(j))/6.0
end do
end subroutine rk4gn

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
      call dPES(flag, a, b, c, d, eo, r, dh, de, n)
      write(1,*) r, h(1), h(2), h(3), h(4), e(1), e(2)
      write(2,*) r, dh(1), dh(2), dh(3), dh(4)
      do while(r<4)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
         call dPES(flag, a, b, c, d, eo, r, dh, de, n)
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
      call dPES(flag, a, b, c, d, eo, r, dh, de, n)
      write(1,*) r, h(1), h(2), h(3), h(4), e(1), e(2)
      write(2,*) r, dh(1), dh(2), dh(3), dh(4)
      do while(r<8)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
         call dPES(flag, a, b, c, d, eo, r, dh, de, n)
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
      call dPES(flag, a, b, c, d, eo, r, dh, de, n)
      write(1,*) r, h(1)*50, h(2), h(3), h(4)*50, e(1), e(2)
      write(2,*) r, dh(1), dh(2), dh(3), dh(4), de(1), de(2)
      do while(r<10)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, n)
         call dPES(flag, a, b, c, d, eo, r, dh, de, n)
         write(1,*) r, h(1)*50, h(2), h(3), h(4)*50, e(1), e(2)
         write(2,*) r, dh(1), dh(2), dh(3), dh(4), de(1), de(2)
      end do
   end if
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt') 
102 format(5(1pe12.3))
end subroutine print_PES

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
