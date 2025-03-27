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
    !call dc_data_collection
    !call ec_data_collection
    !call sb_data_collection
    !call print_pes(n)
    !call MPI_FINALIZE()
end program main

subroutine sc_data_collection
integer, parameter :: n = 4
double precision xi(n+2), avgnumt(n/2), avgnumr(n/2), denom, ip, start, finish
double precision days, hours, minutes, seconds
CALL CPU_TIME(start)
open(unit=1, file='sc_prob.dat')
write(1,*) "# h = 1/(5.0*ip), mc_steps = 15000, gamma = 0.366"
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
write(1,*) '# days = ', days, 'hrs = ', hours, 'min = ', minutes, 's = ', seconds
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
integer reflection ! Flags 1 = yes, 0 = no.
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
gamma      = 0.366!(sqrt(3.0)-1.0)/2.0
dsn        = 2.0*gamma
arg        = 0.0
j          = 0
statef     = 0
numdum     = 1.0
do concurrent (k = 1: n/2)
   numt(k) = 0.0
   numr(k) = 0.0
   wi(k)   = 0.0
   wf(k)   = 0.0
end do

! Parameters.
ti      = 0.0
dt      = 1.0/(5.0*ip)
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
double precision c1, c2, c3
a  = 0.01
b  = 1.6
c  = 0.005
d  = 1.0
mu = 2000.0

c1 = c*exp(-d*x(1)**2.0)
dx(1) = x(2)/mu    ! Rdot
if(x(1) .ge. 0.0) then
   c2 = a*exp(-b*x(1))
   c3 = a - c2
   dx(2) = 2.0*d*x(1)*c1*(x(4)*x(6)+x(3)*x(5)) &
        -c2*b*0.5*(x(3)**2.0+x(4)**2.0-x(5)**2.0-x(6)**2.0)   ! Pdot
   dx(3) = c3*x(4)+x(6)*c1    !x1dot
   dx(4) = -c3*x(3)-x(5)*c1   !p1dot
   dx(5) = x(4)*c1-c3*x(6)    !x2dot
   dx(6) = c3*x(5)-x(3)*c1    !x1dot
else if(x(1) .lt. 0.0) then
   c2 = a*exp(b*x(1))
   c3 = c2 - a
   dx(2) = 2.0*d*x(1)*c1*(x(4)*x(6)+x(3)*x(5)) &
        -c2*b*0.5*(x(3)**2.0+x(4)**2.0-x(5)**2.0-x(6)**2.0)   ! Pdot
   dx(3) = c3*x(4)+x(6)*c1    !x1dot
   dx(4) = -c3*x(3)-x(5)*c1   !p1dot
   dx(5) = x(4)*c1-c3*x(6)    !x2dot
   dx(6) = c3*x(5)-x(3)*c1    !x1dot
end if

!c1 = c*exp(-d*x(1)**2.0) 
!dx(1) = x(2)/mu
!if(x(1) .ge. 0.0) then
!   c2 = a - a*exp(-b*x(1))
!   dx(2) = 2.0*c1*d*x(1)*(x(4)*x(6)+x(3)*x(5))+0.5*(-x(4)**2.0+x(6)**2.0-x(3)**2.0+x(5)**2.0)*a*b*exp(-b*x(1))
!   dx(3) = c1*x(6) + x(4)*c2
!   dx(4) = -c1*x(5) - x(3)*c2
!   dx(5) = c1*x(4) - x(6)*c2
!   dx(6) = -c1*x(3) + x(5)*c2
!else if(x(1) .lt. 0.0) then
!   c2 = a*exp(b*x(1))-a
!   dx(2) = 2.0*c1*d*x(1)*(x(4)*x(6)+x(3)*x(5))+0.5*(-x(4)**2.0+x(6)**2.0-x(3)**2.0+x(5)**2.0)*a*b*exp(b*x(1))
!   dx(3) = c1*x(6) + x(4)*c2
!   dx(4) = -c1*x(5) - x(3)*c2
!   dx(5) = c1*x(4) - x(6)*c2
!   dx(6) = -c1*x(3) + x(5)*c2
!end if
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
write(1,*) "# h = 1/(5.0*ip), mc_steps = 15000, gamma = 0.366"
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
write(1,*) '# days = ', days, 'hrs = ', hours, 'min = ', minutes, 's = ', seconds
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
maxrun = 15000
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
integer reflection ! Flags 1 = yes, 0 = no.
integer statei, statef, k ! E state related.
double precision ti, tf, dt, dumrmin, rmin, rmax, gamma, dsn ! Adjustment parameters.
double precision kdelta, bni(n/2), bnf(n/2), sni(n/2), snf(n/2), q(n/2), rnd(n/2) ! E state related
double precision numt(n/2), numr(n/2), heaviside, arg, wi(n/2), wf(n/2) ! MC Average.
double precision xi(n+2), xf(n+2), ip, energy ! Variables.
double precision numdum
external dc_eq ! Single crossing Hamilton's eqs.
100 continue ! We return after an unsuccessful attempt.
! Initialising variables.
reflection = 0
gamma      = 0.366!(sqrt(3.0)-1.0)/2.0
dsn        = 2.0*gamma
arg        = 0.0
j          = 0
statef     = 0
numdum     = 1.0
do concurrent (k = 1: n/2)
   numt(k) = 0.0
   numr(k) = 0.0
   wi(k)   = 0.0
   wf(k)   = 0.0
end do

! Parameters.
ti      = 0.0
dt      = 1.0/(5.0*ip)
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
double precision c1, c2, c3
a  = 0.1
b  = 0.28
c  = 0.015
d  = 0.06
eo = 0.05
mu = 2000.0
c1 = c*exp(-d*x(1)**2.0)
c2 = a*exp(-b*x(1)**2.0)
c3 = 0.5*(c2-eo)

dx(1) = x(2)/mu    ! Rdot
dx(2) = x(1)*(2.0*d*c1*(x(4)*x(6)+x(3)*x(5))&
        +b*c2*(0.5*(x(3)**2.0+x(4)**2.0-x(5)**2.0-x(6)**2.0)-1)) ! Pdot
dx(3) = c3*x(4)+c1*x(6)  !x1dot
dx(4) = -c3*x(3)-c1*x(5) !p1dot
dx(5) = c1*x(4)-c3*x(6)  !x2dot
dx(6) = c3*x(5)-c1*x(3)  !p2dot

!c1 = 0.5*(a*exp(-b*x(1)**2.0)-eo)
!c2 = c*exp(-d*x(1)**2.0)

!dx(1) = x(2)/mu
!dx(2) = 0.5*exp(-(b+d)*x(1)**2.0)*x(1)*(4*c*d*exp(b*x(1)**2.0)*(x(4)*x(6)+x(3)*x(5))&
!        +a*b*exp(d*x(1)**2.0)*(-2.0+x(4)**2.0-x(6)**2.0+x(3)**2.0-x(5)**2.0))
!dx(3) = c1*x(4)+c2*x(6)
!dx(4) = -c1*x(3)-c2*x(5)
!dx(5) = c2*x(4)-c1*x(6)
!dx(6) = -c2*x(3)+c1*x(5)

!call PES(2, a, b, c, d, eo, x(1), h, e, omega, n)
!call dPES(2, a, b, c, d, eo, x(1), dh, de, n)
!dx(1) = x(2)/mu
!dx(2) = -((dh(1)+dh(4))/2.0 + (x(3)**2.0+x(4)**2.0-x(5)**2.0-x(6)**2.0)*(dh(1)-dh(4))/4.0 + dh(3)*(x(3)*x(5) + x(4)*x(6)))
!dx(3) = x(4)*(h(1)-h(4))/2.0 + x(6)*h(3)
!dx(4) = -(x(3)*(h(1)-h(4))/2.0 + x(5)*h(3))
!dx(5) = -x(6)*(h(1)-h(4))/2.0 + x(4)*h(3)
!dx(6) = -(-x(5)*(h(1)-h(4))/2.0 + x(3)*h(3))
end subroutine dc_eq

subroutine ec_data_collection
integer, parameter :: n = 4
double precision xi(n+2), avgnumt(n/2), avgnumr(n/2), denom, ip, start, finish
double precision days, hours, minutes, seconds
CALL CPU_TIME(start)
open(unit=1, file='ec_prob.dat')
write(1,*) "# h = 1/(5.0*ip), mc_steps = 15000, gamma = 0.366"
write(1,*) "# Initial Momentum  ", "R 1<-1  ", "R 2<-1  ", "T 1<-1  ", "T 2<-1  "

ip = 2.0

!Initial ones.
call ec_mc_avg(avgnumr,avgnumt,denom,ip,n)
write(1,*) ip, avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom

do while(ip < 30.0)
   ip = ip + 1.0
   call ec_mc_avg(avgnumr,avgnumt,denom,ip,n)
   write(1,*) ip, avgnumr(1)/denom, avgnumr(2)/denom, avgnumt(1)/denom, avgnumt(2)/denom
end do

CALL CPU_TIME(finish)
seconds = finish-start
days = floor(seconds/86400)
hours = floor(seconds/3600.0)
minutes = floor(mod(seconds/60.0,60.0))
seconds = mod(seconds,60.0)
write(1,*) '# Execution Time:'
write(1,*) '# days = ', days, 'hrs = ', hours, 'min = ', minutes, 's = ', seconds
end subroutine ec_data_collection

subroutine ec_mc_avg(avgnumr,avgnumt,denom,ip,n)
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
   call ec(ip, numt, numr, n)
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
end subroutine ec_mc_avg

subroutine ec(ip, numt, numr, n)
implicit none
double precision, parameter :: pi = 4.0*atan(1.0)
integer n ! F x F matrix.
integer i, j ! Counters.
integer reflection ! Flags 1 = yes, 0 = no.
integer statei, statef, k ! E state related.
double precision ti, tf, dt, dumrmin, rmin, rmax, gamma, dsn ! Adjustment parameters.
double precision kdelta, bni(n/2), bnf(n/2), sni(n/2), snf(n/2), q(n/2), rnd(n/2) ! E state related
double precision numt(n/2), numr(n/2), heaviside, arg, wi(n/2), wf(n/2) ! MC Average.
double precision xi(n+2), xf(n+2), ip ! Variables.
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
do concurrent (k = 1: n/2)
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
   call rk4gn(ec_eq,ti,tf,xi,xf,n+2)
   ! Check if the particle was reflected back.
   dumrmin = xf(1)
!print*, xi(1)-xf(1)
   !if(xf(1)-xi(1) .lt. 0.0) then
   if(dumrmin .lt. rmin) then
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
double precision t, x(n+2), dx(n+2), h(n), dh(n), e(n), de(n), a, b, c, d, eo, omega, mu
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

subroutine sb_data_collection
integer, parameter :: n = 100
integer j
double precision xi(2*n+4), avgnum(2,21), denom(21), start, finish
double precision days, hours, minutes, seconds
CALL CPU_TIME(start)
open(unit=1, file='sb_prob.dat')
write(1,*) "# h = 1/(5.0*ip), mc_steps = 15000, gamma = 0.366"
write(1,*) "# Initial Momentum  ", "R 1<-1  ", "R 2<-1  ", "T 1<-1  ", "T 2<-1  "

call sb_mc_avg(avgnum,denom,n)

do j = 1, 21
   write(1,*) 0.5*(j-1), (avgnum(1,j)-avgnum(2,j))/denom(j)
end do

CALL CPU_TIME(finish)
seconds = finish-start
days = floor(seconds/86400)
hours = floor(seconds/3600.0)
minutes = floor(mod(seconds/60.0,60.0))
seconds = mod(seconds,60.0)
write(1,*) '# Execution Time:'
write(1,*) '# days = ', days, 'hrs = ', hours, 'min = ', minutes, 's = ', seconds
end subroutine sb_data_collection

subroutine sb_mc_avg(avgnum,denom,n)
implicit none
integer n, i, j, run, maxrun
double precision num(2,21), avgnum(2,21), denom(21)
double precision xi(2*n+4)
! Initialising variables
do concurrent (j = 1: 21)
   do concurrent (i = 1: 2)
      avgnum(i,j) = 0.0
   end do
   denom(j)   = 0.0
end do

! Parameters
maxrun = 1
do concurrent (run = 1: maxrun)
   call sb(num, n)
   ! Summing numerators. Only those corresponding to the final state will change.
   do concurrent (j = 1: 21)
      do concurrent (i = 1: 2)
         avgnum(i,j) = avgnum(i,j) + num(i,j)
      end do
   end do
end do

do concurrent (j = 1: 21)
   do concurrent (i = 1: 2)
      avgnum(i,j) = avgnum(i,j)/maxrun
      denom(j) = denom(j) + avgnum(i,j)
   end do
end do
end subroutine sb_mc_avg

subroutine sb(num, n)
implicit none
double precision, parameter :: pi = 4.0*atan(1.0)
integer n, k, j, i, statei, statef
double precision xi(2*n+4), xf(2*n+4), w(n), c(n), dw, wc, m, a, b, gamma, rnd, rnd2(2)
double precision bni(2), sni(2), wi(2), bnf(2), snf(2), wf(2), q(2), kdelta, heaviside, arg, dsn
double precision num(2,21), numdum, ti, dt, tf, tmax
double precision c1, c2, c3
! Normal distribution RNG data declarations
  real ( kind = 4 ) fn(128)
  integer ( kind = 4 ) kn(128)
  real ( kind = 4 ) r4_nor
  integer ( kind = 4 ) seed
  real ( kind = 8 ) rnrm
  real ( kind = 4 ) wn(128)
  double precision r
!------------------------------------------
external sb_eq ! Equations of motion.
! General variables
100 continue
gamma      = 0.366!(sqrt(3.0)-1.0)/2.0
dsn        = 2.0*gamma
arg        = 0.0
statef     = 0
numdum     = 1.0
do concurrent (k = 1: 2)
   num(k,1)  = 0.0 ! Initial state
   wi(k)   = 0.0
   wf(k)   = 0.0
end do

! Parameters.
ti      = 0.0
dt      = 0.01
tmax    = 10.0
statei  = 1

! Spin boson variables
wc       = 1.0     ! Characteristic frequency.
w(:)     = 0.0     ! Oscillating frequencies.
m        = 2000.0  ! Particle mass.
a        = 0.09    ! Kondo parameter.
b        = 12.5

! Assigning oscillating frequencies.
j = 1 ! Counter.
call init_random_seed()
do while(j<=n)
   call random_number(rnd)
   w(j) = rnd*4.0*wc
   if(w(j) .ge. 0.01*wc .and. w(j) .le. 4.0*wc) then
      j = j + 1
   end if
end do
j = 0

! Initial electronic states.
do k = 1, 2
   call random_number(rnd2)
   bni(k) = kdelta(statei,k)
   sni(k) = bni(k) + gamma*(2.0*rnd2(1)-1.0)
   q(k)   = 2.0*pi*rnd2(2)
   arg    = gamma-abs(sni(k)-bni(k))
   wi(k)  = heaviside(arg)/dsn
end do

dw = (maxval(w)-minval(w))/(n-1.0) ! Delta_w

! Calculating proportionality coefficients.
do concurrent (k = 1:n)
   c(k) = w(k)*sqrt(a*dw*m*exp(w(k)/wc))
end do

! Next comment section is no longer relevant, the problem has been solved.
!---------------------------------------------------------------------!
! For initial parameter sampling from the wigner distribution use:    !
! http://www.gnu.org/software/mcsim/                                  !
! We have to load the initial parameters from the text file generated !
! [Load text file and assign to data(:,:)]                            !
!---------------------------------------------------------------------!

call r4_nor_setup ( kn, fn, wn )
call random_number(r)
seed = floor(r*1000000000)
do k = 1,n
   c1 = sqrt(2.0*tanh(b*w(k)/2.0))
   c2 = w(k)*m
   c3 = sqrt(c2)
   ! Initial nuclear positions
   rnrm = r4_nor ( seed, kn, fn, wn )
   xi(k) = rnrm*c3/c1  
   ! Initial nuclear momenta
   rnrm = r4_nor ( seed, kn, fn, wn )
   xi(k+n) = rnrm/(c1*c3) - c(k)/(c2*w(k))
end do

! Initial Positions.
!do concurrent (k = 1:100)
!   xi(k) = data(k,1)
!end do
! Initial Momenta.
!do concurrent (k = 101:200)
!   xi(k) = data(k,2)
!end do
! Initial electronic coordinates
xi(201) = sqrt(2.0*(sni(1)+gamma))*cos(q(1))  ! x1
xi(202) = -sqrt(2.0*(sni(1)+gamma))*sin(q(1)) ! p1
xi(203) = sqrt(2.0*(sni(2)+gamma))*cos(q(2))  ! x2
xi(204) = -sqrt(2.0*(sni(2)+gamma))*sin(q(2)) ! p2

! Numerator for the initial population states.
j = 1
do concurrent (k = 1: 2)
      numdum = numdum*wi(k)*wi(k)
end do
num(statef,j) = numdum
numdum = 1.0 ! Reset for the next time it's used.

! Integrate hamilton's equations
do while(ti<=tmax)
   tf = ti + dt
   call rk4gn(sb_eq,ti,tf,xi,xf,2*n+4)
   ! Prepares next iteration
   ti = tf
   do concurrent (i = 1: 2*n+4)
      xi(i) = xf(i)
   end do
   ! Calculating intermediate populations.
   if (mod(tf,0.5) .eq. 0.0) then ! Every 0.5 units of time.
      j = j + 1 ! Counter for intermediate time step.
      i = 0     ! counter used to calculate nk(t).
      ! Calculating nk(t).
      do concurrent (k = 1: 2)
         snf(k) = 0.5*xf(2+k+i)**2.0 + 0.5*xf(2+k+1+i)**2.0
         i      = i + 1
      end do
      if(1.0 .le. snf(1) .and. snf(1) .le. 1.0 + dsn .and. 0.0 .le. snf(2) .and. snf(2) .le. dsn) then
         statef  = 1
      else if(0.0 .le. snf(1) .and. snf(1) .le. dsn .and. 1.0 .le. snf(2) .and. snf(2) .le. 1.0 + dsn) then
         statef  = 2
      else
         ! Well, go to 100 (the start of the subroutine).
         go to 100
      end if
      ! Intermediate Electronic states.
      do k = 1, 2
         bnf(k) = kdelta(statef,k)
         arg    = gamma - abs(snf(k)-bnf(k))
         wf(k)  = heaviside(arg)/dsn
      end do
      do concurrent (k = 1: 2)
         numdum = numdum*wi(k)*wf(k)
      end do
      ! Only the numerator which corresponds to the final state will change, the other will remain equal to 0.0
      num(statef,j) = numdum
      numdum = 1.0
   end if
end do
end subroutine sb

subroutine sb_eq(t,x,dx,w,c,n)
!===========================================================================!
! Equations of motion of single avoided crossing                            !
! Daniel Celis Garza 25 March 2015                                          !
!---------------------------------------------------------------------------!
! x(1)   -> x(100) = R1 -> R100     dx(1)   -> dx(100) = dR1/dt -> dR100/dt !
! x(101) -> x(200) = P1 -> P100     dx(101) -> dx(200) = dP1/dt -> dP100/dt !
! x(201) = x1                       dx(3)              = dx1/dt             !
! x(202) = p1                       dx(4)              = dp1/dt             !
! x(203) = x2                       dx(5)              = dx2/dt             !
! x(204) = p2                       dx(6)              = dp2/dt             !
!---------------------------------------------------------------------------!
! h(1) = H11  h(2) = H12  h(3) = H21  h(4) = H22                            !
! H11, H22, H33, H44 have to be differentiated individually.                !
!===========================================================================!
integer n, i
double precision t, x(2*n+4), dx(2*n+4), w(n), c(n), m, v0, v1, eps, delta
double precision c1, dc1
m  = 2000.0
eps = 1.0
delta = 1.0/2.5
v0 = 0.0
v1 = 0.0

! V0 and V1
do concurrent (i = 1:100)
   v0 = v0 + w(i)**2*x(i)
   v1 = v1 + c(i)*x(i)
end do
v0 = 1000.0*v0

! H11, H22, H33 and H44
!h(1) = v0 + v1 + eps
!h(2) = delta
!h(3) = delta
!h(4) = v0 - v1 - eps

c1 = v1 + eps
! Equations of motion.
! Positions and momenta.
do concurrent (i = 1:100)
   dx(i)     = x(i+100)/m
   dx(i+100) = -m*w(i)**2.0*x(i) &
               - 0.5*c(i)*(x(201)**2.0+x(202)**2.0-x(203)**2.0-x(204)**2.0)
end do
dx(201) = x(202)*c1 + x(204)*delta
dx(202) = -x(201)*c1 - x(203)*delta
dx(203) = -x(204)*c1 + x(202)*delta
dx(204) = x(203)*c1 - x(201)*delta
end subroutine sb_eq

! real function rnd_exp()
!implicit none
!real rnd
!call random_number(rnd)
!rnd_exp = -log(rnd)
!return
!end function rnd_exp

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
   e(1) = exp(-(b+d)*r**2)*(-a*exp(d*r**2)+exp((b+d)*r**2)*eo-sqrt(4*c**2*exp(2*b*r**2)&
        +exp(2*d*r**2)*(a-exp(b*r**2)*eo)**2))/2
   e(2) = exp(-(b+d)*r**2)*(-a*exp(d*r**2)+exp((b+d)*r**2)*eo+sqrt(4*c**2*exp(2*b*r**2)&
        +exp(2*d*r**2)*(a-exp(b*r**2)*eo)**2))/2
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

!****************************************************
! Ziggurat algorithm
! I modified it so that the seeds are random not fixed.
!****************************************************
function r4_exp ( jsr, ke, fe, we )

!*****************************************************************************80
!
!! R4_EXP returns an exponentially distributed single precision real value.
!
!  Discussion:
!
!    The underlying algorithm is the ziggurat method.
!
!    Before the first call to this function, the user must call R4_EXP_SETUP
!    to determine the values of KE, FE and WE.
!
!    Note that JSR and KE should actually be of unsigned integer type, but
!    this data type is not supported by FORTRAN90.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Input, integer ( kind = 4 ) KE(256), data computed by R4_EXP_SETUP.
!
!    Input, real ( kind = 4 ) FE(256), WE(256), data computed by R4_EXP_SETUP.
!
!    Output, real ( kind = 4 ) R4_EXP, an exponentially distributed
!    random value.
!
  implicit none

  real ( kind = 4 ) fe(256)
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) ke(256)
  real ( kind = 4 ) r4_exp
  real ( kind = 4 ) r4_uni
  integer ( kind = 4 ) shr3
  real ( kind = 4 ) value
  real ( kind = 4 ) we(256)
  real ( kind = 4 ) x

  jz = shr3 ( jsr )
  iz = iand ( jz, 255 )

  if ( abs ( jz  ) < ke(iz+1) ) then

    value = real ( abs ( jz ), kind = 4 ) * we(iz+1)

  else

    do

      if ( iz == 0 ) then
        value = 7.69711E+00 - log ( r4_uni ( jsr ) )
        exit
      end if

      x = real ( abs ( jz ), kind = 4 ) * we(iz+1)

      if ( fe(iz+1) + r4_uni ( jsr ) * ( fe(iz) - fe(iz+1) ) &
           < exp ( - x ) )  then
        value = x
        exit
      end if

      jz = shr3 ( jsr )
      iz = iand ( jz, 255 )

      if ( abs ( jz ) < ke(iz+1) ) then
        value = real ( abs ( jz ), kind = 4 ) * we(iz+1)
        exit
      end if

    end do

  end if

  r4_exp = value

  return
end
subroutine r4_exp_setup ( ke, fe, we )

!*****************************************************************************80
!
!! R4_EXP_SETUP sets data needed by R4_EXP.
!
!  Discussion:
!
!    Note that KE should actually be of unsigned integer type, but
!    this data type is not supported by FORTRAN90.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KE(256), data needed by R4_EXP.
!
!    Output, real ( kind = 4 ) FE(256), WE(256), data needed by R4_EXP.
!
  implicit none

  real ( kind = 8 ) de
  real ( kind = 4 ) fe(256)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ke(256)
! real ( kind = 8 ), parameter :: m2 = 4294967296.0D+00
  real ( kind = 8 ), parameter :: m2 = 2147483648.0D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) te
  real ( kind = 8 ), parameter :: ve = 3.949659822581572D-03
  real ( kind = 4 ) we(256)

  de = 7.697117470131487D+00
  te = 7.697117470131487D+00

  q = ve / exp ( - de )

  ke(1) = int ( ( de / q ) * m2 )
  ke(2) = 0

  we(1) = real ( q / m2, kind = 4 )
  we(256) = real ( de / m2, kind = 4 )

  fe(1) = 1.0E+00
  fe(256) = real ( exp ( - de ), kind = 4 )

  do i = 255, 2, -1
    de = - log ( ve / de + exp ( - de ) )
    ke(i+1) = int ( ( de / te ) * m2 )
    te = de
    fe(i) = real ( exp ( - de ), kind = 4 )
    we(i) = real ( de / m2, kind = 4 )
  end do

  return
end
function r4_nor ( jsr, kn, fn, wn )

!*****************************************************************************80
!
!! R4_NOR returns a normally distributed single precision real value.
!
!  Discussion:
!
!    The value returned is generated from a distribution with mean 0 and
!    variance 1.
!
!    The underlying algorithm is the ziggurat method.
!
!    Before the first call to this function, the user must call R4_NOR_SETUP
!    to determine the values of KN, FN and WN.
!
!    Note that JSR and KN should actually be of unsigned integer type, but
!    this data type is not supported by FORTRAN90.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Input, integer ( kind = 4 ) KN(128), data computed by R4_NOR_SETUP.
!
!    Input, real ( kind = 4 ) FN(128), WN(128), data computed by R4_NOR_SETUP.
!
!    Output, real ( kind = 4 ) R4_NOR, a normally distributed random value.
!
  implicit none

  real ( kind = 4 ) fn(128)
  integer ( kind = 4 ) hz
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) kn(128)
  real ( kind = 4 ), parameter :: r = 3.442620E+00
  real ( kind = 4 ) r4_nor
  real ( kind = 4 ) r4_uni
  integer ( kind = 4 ) shr3
  real ( kind = 4 ) value
  real ( kind = 4 ) wn(128)
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  hz = shr3 ( jsr )
  iz = iand ( hz, 127 )

  if ( abs ( hz ) < kn(iz+1) ) then

    value = real ( hz, kind = 4 ) * wn(iz+1)

  else

    do

      if ( iz == 0 ) then

        do
          x = - 0.2904764E+00 * log ( r4_uni ( jsr ) )
          y = - log ( r4_uni ( jsr ) )
          if ( x * x <= y + y ) then
            exit
          end if
        end do

        if ( hz <= 0 ) then
          value = - r - x
        else
          value = + r + x
        end if

        exit

      end if

      x = real ( hz, kind = 4 ) * wn(iz+1)

      if ( fn(iz+1) + r4_uni ( jsr ) * ( fn(iz) - fn(iz+1) ) &
         < exp ( - 0.5E+00 * x * x ) ) then
        value = x
        exit
      end if

      hz = shr3 ( jsr )
      iz = iand ( hz, 127 )

      if ( abs ( hz ) < kn(iz+1) ) then
        value = real ( hz, kind = 4 ) * wn(iz+1)
        exit
      end if

    end do

  end if

  r4_nor = value

  return
end
subroutine r4_nor_setup ( kn, fn, wn )

!*****************************************************************************80
!
!! R4_NOR_SETUP sets data needed by R4_NOR.
!
!  Discussion:
!
!    Note that KN should actually be of unsigned integer type, but
!    this data type is not supported by FORTRAN90.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KN(128), data needed by R4_NOR.
!
!    Output, real ( kind = 4 ) FN(128), WN(128), data needed by R4_NOR.
!
  implicit none

  real ( kind = 8 ) dn
  real ( kind = 4 ) fn(128)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kn(128)
  real ( kind = 8 ), parameter :: m1 = 2147483648.0D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) tn
  real ( kind = 8 ), parameter :: vn = 9.91256303526217D-03
  real ( kind = 4 ) wn(128)

  dn = 3.442619855899D+00
  tn = 3.442619855899D+00

  q = vn / exp ( - 0.5D+00 * dn * dn )

  kn(1) = int ( ( dn / q ) * m1 )
  kn(2) = 0

  wn(1) = real ( q / m1, kind = 4 )
  wn(128) = real ( dn / m1, kind = 4 )

  fn(1) = 1.0E+00
  fn(128) = real ( exp ( - 0.5D+00 * dn * dn ), kind = 4 )

  do i = 127, 2, -1
    dn = sqrt ( - 2.0D+00 * log ( vn / dn + exp ( - 0.5D+00 * dn * dn ) ) )
    kn(i+1) = int ( ( dn / tn ) * m1 )
    tn = dn
    fn(i) = real ( exp ( - 0.5D+00 * dn * dn ), kind = 4 )
    wn(i) = real ( dn / m1, kind = 4 )
  end do

  return
end
function r4_uni ( jsr )

!*****************************************************************************80
!
!! R4_UNI returns a uniformly distributed real value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Output, real ( kind = 4 ) R4_UNI, a uniformly distributed random value in
!    the range [0,1].
!
  implicit none

  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_input
  real ( kind = 4 ) r4_uni

  jsr_input = jsr

  jsr = ieor ( jsr, ishft ( jsr,   13 ) )
  jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
  jsr = ieor ( jsr, ishft ( jsr,    5 ) )

! r4_uni = 0.5E+00 + 0.2328306E-09 * real ( jsr_input + jsr, kind = 4 )

  r4_uni = 0.5E+00 + real ( jsr_input + jsr, kind = 4 ) &
    / real ( 65536, kind = 4 ) / real ( 65536, kind = 4 )

  return
end
function shr3 ( jsr )

!*****************************************************************************80
!
!! SHR3 evaluates the SHR3 generator for integers.
!
!  Discussion:
!
!    Note that JSR should actually be of unsigned integer type, but
!    this data type is not supported by FORTRAN77.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Output, integer ( kind = 4 ) SHR3, the value of the SHR3 generator.
!
  implicit none

  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_input
  integer ( kind = 4 ) shr3

  jsr_input = jsr

  jsr = ieor ( jsr, ishft ( jsr,   13 ) )
  jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
  jsr = ieor ( jsr, ishft ( jsr,    5 ) )

  shr3 = jsr_input + jsr

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

!**********************************
! Random seed
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
