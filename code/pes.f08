program main
  implicit none
  integer, parameter :: es = 2 ! Electronic states.
  call print_pes(es)
end program main

subroutine PES(flag, a, b, c, d, eo, r, h, e, omega, es)
!=========================================================!
! Potential energy surfaces as defined by Miller          !
!---------------------------------------------------------!
! Daniel Celis Garza 16 Jan 2015, 20 Jan 2015             !
!---------------------------------------------------------!
! flag = 1  ->  single avoided crossing (diabatic)        !
! flag = 2  ->  dual avoided crossing (diabatic)          !
! flag = 3  ->  extended coupling (adiabatic)             !
!---------------------------------------------------------!
! h(n)      -> diabatic matrix elements                   !
! e(n/2)    -> adiabatic matrix elements (eigenvalues of  !
!              the diabatic matrix)                       !
!---------------------------------------------------------!
! h(1) = H11    h(2) = H12    h(3) = H21    h(4) = H22    !
! e(1) = E1     e(2) = E2                                 !
!=========================================================!
implicit none
integer es
integer flag
double precision a, b, c, d, eo, r, h(es*es), e(es), omega
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

subroutine dPES(flag, a, b, c, d, r, dh, de, es)
!===============================================================!
! Analytical derivatives of the PES used in Hamilton's eqs      !
!---------------------------------------------------------------!
! Daniel Celis Garza 21 Jan 2015                                !
!===============================================================!
implicit none
integer es
integer flag
double precision a, b, c, d, r, dh(es*es), de(es)
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

subroutine print_PES(es)
! Testing testing all PES
implicit none
integer es
integer flag
double precision a, b, c, d, eo, r, h(es*es), e(es), omega, dh(es*es), de(es)
do flag = 1, 3
   if (flag == 1) then
      r = -4.0
      a = 0.01
      b = 1.6
      c = 0.005
      d = 1
      open(unit=1, file='sac.dat')
      open(unit=2, file='dsac.dat')
      call PES(flag, a, b, c, d, eo, r, h, e, omega, es)
      call dPES(flag, a, b, c, d, r, dh, de, es)
      write(1,*) r, h(1), h(2), h(3), h(4), e(1), e(2)
      write(2,*) r, dh(1), dh(2), dh(3), dh(4)
      do while(r<4)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, es)
         call dPES(flag, a, b, c, d, r, dh, de, es)
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
      call PES(flag, a, b, c, d, eo, r, h, e, omega, es)
      call dPES(flag, a, b, c, d, r, dh, de, es)
      write(1,*) r, h(1), h(2), h(3), h(4), e(1), e(2)
      write(2,*) r, dh(1), dh(2), dh(3), dh(4)
      do while(r<8)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, es)
         call dPES(flag, a, b, c, d, r, dh, de, es)
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
      call PES(flag, a, b, c, d, eo, r, h, e, omega, es)
      call dPES(flag, a, b, c, d, r, dh, de, es)
      write(1,*) r, h(1)*50, h(2), h(3), h(4)*50, e(1), e(2)
      write(2,*) r, dh(1), dh(2), dh(3), dh(4), de(1), de(2)
      do while(r<10)
         r = r + 0.01
         call PES(flag, a, b, c, d, eo, r, h, e, omega, es)
         call dPES(flag, a, b, c, d, r, dh, de, es)
         write(1,*) r, h(1)*50, h(2), h(3), h(4)*50, e(1), e(2)
         write(2,*) r, dh(1), dh(2), dh(3), dh(4), de(1), de(2)
      end do
   end if
end do
end subroutine print_PES
