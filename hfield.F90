subroutine hfield()
    use constants
    implicit none
    integer :: i
  
    !$omp parallel do
    do i=istart,iend
        hx(i) = amx(i)*hx(i)
    enddo
    !$omp end parallel do

    call exchg1d(ez,istart,iend,comm,left,right)
    !$omp parallel do
    do i=istart,iend
        hy(i) = amy(i)*hy(i)+bmyx(i)*(ez(i+1)-ez(i))
    enddo
    !$omp end parallel do

    call exchg1d(ey,istart,iend,comm,left,right)
    !$omp parallel do
    do i=istart,iend
        hz(i) = amz(i)*hz(i)-bmzx(i)*(ey(i+1)-ey(i))
    enddo
    !$omp end parallel do

end subroutine