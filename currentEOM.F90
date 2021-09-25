subroutine currentEOM()
    use constants
    use omp_lib
    implicit none
    integer :: i

    !$omp parallel do
    do i=0,nx
        jx(i) = -qe*nd(i)*vx(i)
        jy(i) = -qe*nd(i)*vy(i)
        jz(i) = -qe*nd(i)*vz(i)
    enddo
    !$omp end parallel do
            
end subroutine 