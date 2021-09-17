subroutine phi()
    use constants
    implicit none
    integer :: i

    !$omp parallel do
    do i=0,nx-1
        phix(i) = (dxhi0-dxi0)*ex(i)+dxi0*pex(i)+exp(-nu*dt)*phix(i) 
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,nx-1
        phiy(i) = (dxhi0-dxi0)*ey(i)+dxi0*pey(i)+exp(-nu*dt)*phiy(i) 
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,nx-1
        phiz(i) = (dxhi0-dxi0)*ez(i)+dxi0*pez(i)+exp(-nu*dt)*phiz(i) 
    enddo
    !$omp end parallel do
     
end subroutine