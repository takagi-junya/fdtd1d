subroutine hfield()
    use constants
    implicit none
    integer :: i
  
    !$omp parallel do
    do i=1,nx-1
        hx(i) = amx(i)*hx(i)-bmxy(i)*(ez(i)-ez(i))
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=0,nx-1
        hy(i) = amy(i)*hy(i)+bmyx(i)*(ez(i+1)-ez(i))
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=0,nx-1
        hz(i) = amz(i)*hz(i)&
        &       -bmzx(i)*(ey(i+1)-ey(i))&
        &       +bmzy(i)*(ex(i)-ex(i))
    enddo
    !$omp end parallel do

end subroutine