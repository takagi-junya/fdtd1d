subroutine efield()
    use constants
    implicit none
    integer :: i

    !$omp parallel do
    do i=istart,iend
        ex(i) = aex(i)*ex(i)
        ex(i) = ex(i) - ajj*jx(i)
    enddo
    !$omp end parallel do

    call exchg1d(hz,istart,iend,comm,left,right)
    !$omp parallel do
    do i=istart,iend
        ey(i) = aey(i)*ey(i)-beyx(i)*(hz(i)-hz(i-1))
        ey(i) = ey(i) - ajj*jy(i)
    enddo
    !$omp end parallel do

    call exchg1d(hy,istart,iend,comm,left,right)
    !$omp parallel do
    do i=istart,iend
        ez(i) = aez(i)*ez(i)+bezx(i)*(hy(i)-hy(i-1))
        ez(i) = ez(i) - ajj*jz(i)
    enddo
    !$omp end parallel do
     
end subroutine