subroutine efield()
    use constants
    implicit none
    integer :: i

    !$omp parallel do
    do i=0,nx-1
        ex(i) = aex(i)*ex(i) + bexy(i)*(hz(i)-hz(i))
        ex(i) = ex(i) - ajj*jx(i)
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,nx-1
        ey(i) = aey(i)*ey(i) - beyx(i)*(hz(i)-hz(i-1))
        ey(i) = ey(i) - ajj*jy(i)
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,nx-1
        ez(i) = aez(i)*ez(i)&
        &       + bezx(i)*(hy(i)-hy(i-1))&
        &       - bezy(i)*(hx(i)-hx(i))
        ez(i) = ez(i) - ajj*jz(i)
    enddo
    !$omp end parallel do
     
end subroutine