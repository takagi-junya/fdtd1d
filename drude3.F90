subroutine drude3()
    use constants 
    implicit none
    integer ::i
    ajj = dt/eps0 
    do i=ic-prad,ic+prad-1
        ajx(i) = (2.0d0-nu*dt)/(2.0d0+nu*dt)
        ajex(i)= (2.0d0*eps0*(wp**2.0d0)*dt)/(2.0d0+nu*dt)
    enddo

    do i=ic-prad,ic+prad
        ajy(i) = (2.0d0-nu*dt)/(2.0d0+nu*dt)
        ajey(i) = (2.0d0*eps0*(wp**2.0d0)*dt)/(2.0d0+nu*dt)
    enddo

    do i=ic-prad,ic+prad
        ajz(i) = (2.0d0-nu*dt)/(2.0d0+nu*dt)
        ajez(i) = (2.0d0*eps0*(wp**2.0d0)*dt)/(2.0d0+nu*dt)
    enddo
    
end subroutine
            