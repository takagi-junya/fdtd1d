subroutine efield_plrc()
    use constants
    implicit none
    integer :: i

    !$omp parallel do
    do i=0,nx-1
        pex(i) = ex(i)
        ex(i) = aex(i)*ex(i) + bexy(i)*(hz(i)-hz(i))
        ex(i) = ex(i) + aphix(i)*phix(i)
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,nx-1
        pey(i) = ey(i)
        ey(i) = aey(i)*ey(i) - beyx(i)*(hz(i)-hz(i-1))
        ey(i) = ey(i) + aphiy(i)*phiy(i)
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,nx-1
        pez(i) = ez(i)
        ez(i) = aez(i)*ez(i)&
        &       + bezx(i)*(hy(i)-hy(i-1))&
        &       - bezy(i)*(hx(i)-hx(i))
        ez(i) = ez(i) + aphiz(i)*phiz(i)
    enddo
    !$omp end parallel do
     
end subroutine