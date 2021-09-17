module pwave
    implicit none
    real :: orge(3),orgb(3)
contains

    subroutine init_pwave
        use constants
        implicit none
        integer :: i
        real :: x
        t=0.0
        do i=1,nx-1
            x = i*dx
            ey(i) = amps(2)*gs_ez(x)
        enddo

        do i=1,nx-1
            x = i*dx
            ez(i) = amps(3)*gs_ez(x)
        enddo

    t=t+0.5*dt
        do i=0,nx-1
            x = (i+0.5)*dx
            hy(i) = -amps(3)*gs_ez(x)/z0
        enddo
        do i=1,nx-1
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
        gs_ez = -2.0d0*a*pw*exp(-a*a)
    end function gs_ez

end module pwave

