subroutine velocity
    use constants
    implicit none
    integer :: i

    if(iple.eq.ic+width) then
        iple = iple-1
    endif

    !$omp parallel do
    do i=ipls,iple
        vx(i) = sab(1,1)*vx(i)+sab(1,2)*vy(i)+sab(1,3)*vz(i)-tc(1,1)*ex(i)-tc(1,2)*ey(i)-tc(1,3)*ez(i)
    enddo
    !$omp end parallel do

    if(iple.eq.ic+width-1) then
        iple = iple+1 
    endif

    !$omp parallel do
    do i=ipls,iple
        vy(i) = sab(2,1)*vx(i)+sab(2,2)*vy(i)+sab(2,3)*vz(i)-tc(2,1)*ex(i)-tc(2,2)*ey(i)-tc(2,3)*ez(i)
        vz(i) = sab(3,1)*vx(i)+sab(3,2)*vy(i)+sab(3,3)*vz(i)-tc(3,1)*ex(i)-tc(3,2)*ey(i)-tc(3,3)*ez(i)
    enddo
    !$omp end parallel do

end subroutine 
            
