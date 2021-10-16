subroutine ADE()
    use constants 
    implicit none
    integer ::i
    ajj = dt/eps0 
    if(iple.eq.ic+width) then
        iple = iple - 1 
    endif
    do i=ipls,iple
        ajx(i) = (2.0d0-nu*dt)/(2.0d0+nu*dt)
        ajex(i)= (2.0d0*eps0*(wp**2.0d0)*dt)/(2.0d0+nu*dt)
    enddo

    if(iple.eq.ic+width-1) then
        iple = iple + 1 
    endif
    do i=ipls,iple
        ajy(i) = (2.0d0-nu*dt)/(2.0d0+nu*dt)
        ajey(i) = (2.0d0*eps0*(wp**2.0d0)*dt)/(2.0d0+nu*dt)
    enddo

    do i=ipls,iple
        ajz(i) = (2.0d0-nu*dt)/(2.0d0+nu*dt)
        ajez(i) = (2.0d0*eps0*(wp**2.0d0)*dt)/(2.0d0+nu*dt)
    enddo
    
end subroutine
            