module functions
implicit none
contains
  function heaviside(x)
    !=======================================================!
    ! Heaviside function                                    !
    ! Daniel Celis Garza 15 Jan 2015                        !
    !=======================================================!
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
    integer i, j
    double precision kdelta
    kdelta = 0.0
    if (i .eq. j) kdelta = 1.0
  end function kdelta

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
    integer n, j
    double precision ti, tf, h, t, xi(n), xf(n), k1(n), k2(n), k3(n), k4(n), x(n), dx(n), sqr2

    h = tf - ti ! This will be replaced by the error-correcting routine.
    t = ti
    sqr2 = sqrt(2.0)
    ! Calculate k1
    call f(t,xi,dx,n)            ! Call the equations of motion
    do concurrent (j = 1: n)     ! Go through all equations of motion
       k1(j) = h*dx(j)           ! Calculate the value of k1 for each equation of motion
       x(j)  = xi(j) + k1(j)/2.0 ! Calculate the next value of x for each equation
    end do                       ! (to be used in the next ki)

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

  subroutine rk4gnsb(f,ti,tf,xi,xf,w,c,n,es)
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
    integer n, es, j
    double precision ti,tf,h,t,xi(n),xf(n),k1(n),k2(n),k3(n),k4(n),x(n),dx(n),sqr2,w(n),c(n)

    h = tf - ti ! This will be replaced by the error-correcting routine.
    t = ti
    sqr2 = sqrt(2.0)
    ! Calculate k1
    call f(t,xi,dx,w,c,n,es)     ! Call the equations of motion
    do concurrent (j = 1: n)     ! Go through all equations of motion
       k1(j) = h*dx(j)           ! Calculate the value of k1 for each equation of motion
       x(j)  = xi(j) + k1(j)/2.0 ! Calculate the next value of x for each equation
    end do                       ! (to be used in the next ki)

    ! Calculate k2
    call f(t+h/2.0,x,dx,w,c,n,es)
    do concurrent (j = 1: n)
       k2(j) = h*dx(j)
       x(j)  = xi(j) + (sqr2-1.0)*k1(j)/2.0 + (1.0-sqr2/2.0)*k2(j)
    end do

    ! Calculate k3
    call f(t+h/2.0,x,dx,w,c,n,es)
    do concurrent (j = 1: n)
       k3(j) = h*dx(j)
       x(j)  = xi(j) - sqr2*k2(j)/2.0 + (1.0+sqr2/2.0)*k3(j)
    end do

    ! Calculate k4 and xf
    call f(t+h,x,dx,w,c,n,es)
    do concurrent (j = 1: n)
       k4(j) = h*dx(j)
       xf(j)  = xi(j) + (k1(j) + (2.0-sqr2)*k2(j) + (2+sqr2)*k3(j) + k4(j))/6.0
    end do
  end subroutine rk4gnsb

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
end module functions
