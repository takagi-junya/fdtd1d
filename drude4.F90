subroutine drude4()
    use constants 
    implicit none
    integer ::i
    ajj = dt/eps0 
    do i=ic-prad,ic+prad-1
        ajx(i)  = exp(-nu*dt)
        ajex(i) = eps0*(wp**2.0d0)*exp(-nu*dt*0.50d0)*dt
    enddo

    do i=ic-prad,ic+prad
        ajy(i)  = exp(-nu*dt)
        ajey(i) = eps0*(wp**2.0d0)*exp(-nu*dt*0.50d0)*dt
    enddo

    do i=ic-prad,ic+prad
        ajz(i)  = exp(-nu*dt)
        ajez(i) = eps0*(wp**2.0d0)*exp(-nu*dt*0.50d0)*dt
    enddo
    
end subroutine
            