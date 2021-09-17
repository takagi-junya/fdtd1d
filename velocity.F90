subroutine velocity
    use constants
    implicit none
    integer :: i

    !$omp parallel do
    do i=ic-prad,ic+prad-1
        !vx(i) = avx(i)*vx(i)+ajex(i)*ex(i)
        vx(i) = omata(1,1)*vx(i)+omata(1,2)*vy(i)+omata(1,3)*vz(i)-omatb(1,1)*ex(i)-omatb(1,2)*ey(i)-omatb(1,3)*ez(i)
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=ic-prad,ic+prad
        vy(i) = omata(2,1)*vx(i)+omata(2,2)*vy(i)+omata(2,3)*vz(i)-omatb(2,1)*ex(i)-omatb(2,2)*ey(i)-omatb(2,3)*ez(i)
        vz(i) = omata(3,1)*vx(i)+omata(3,2)*vy(i)+omata(3,3)*vz(i)-omatb(3,1)*ex(i)-omatb(3,2)*ey(i)-omatb(3,3)*ez(i)
    enddo
    !$omp end parallel do

end subroutine 
            
