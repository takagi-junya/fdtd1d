!-----------------------------------------
! 　　左向き円偏波の設定
!-----------------------------------------
module LCPwave
    implicit none
    real :: orge(3),orgb(3)
contains

    subroutine LCP
        use constants
        implicit none
        integer :: i
        real :: x
        t=0.0
        do i=istart,iend
            x = i*dx
            ey(i) = amps(2)*gs_ez(x)
        enddo

        do i=istart,iend
            x = i*dx
            ez(i) = amps(3)*gs_ez2(x)
        enddo

        t=t+0.5*dt
        do i=istart,iend
            x = (i+0.5)*dx
            hy(i) = -amps(3)*gs_ez2(x)/z0
        enddo
        do i=istart,iend
            x = (i+0.5)*dx
            hz(i) = amps(2)*gs_ez(x)/z0
        enddo

    end subroutine

    real function gs_ez(x)
        use constants
        implicit none
        real,intent(in) :: x
        real :: xx,a
        xx = x-c*t-pc
        a = xx/pw
        gs_ez = sin(omega*t-k0*x)*exp(-a*a)
    end function gs_ez

    real function gs_ez2(x)
        use constants
        implicit none
        real,intent(in) :: x
        real :: xx,a
        xx = x-c*t-pc
        a = xx/pw
        gs_ez2 = cos(omega*t-k0*x)*exp(-a*a)
    end function gs_ez2

end module