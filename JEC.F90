subroutine JEC()
    use constants 
    implicit none
    integer ::i
    ajj = dt/eps0 
    if(iple.eq.ic+width) then
        iple = iple - 1 
    endif 
    do i=ipls,iple
        ajx(i)  = exp(-nu*dt)
        ajex(i) = eps0*(wp**2.0d0)*exp(-nu*dt*0.50d0)*dt
    enddo

    if(iple.eq.ic+width-1) then
        iple = iple + 1 
    endif 

    do i=ipls,iple
        ajy(i)  = exp(-nu*dt)
        ajey(i) = eps0*(wp**2.0d0)*exp(-nu*dt*0.50d0)*dt
    enddo

    do i=ipls,iple
        ajz(i)  = exp(-nu*dt)
        ajez(i) = eps0*(wp**2.0d0)*exp(-nu*dt*0.50d0)*dt
    enddo
    
end subroutine
            