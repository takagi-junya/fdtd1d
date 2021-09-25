!------------------------------------------------------
!   calclate current density caused by plasma
!------------------------------------------------------
subroutine current()
    use constants
    use omp_lib
    implicit none
    integer i
        !$omp parallel do
        do i=0,nx
            jx(i) = ajx(i)*jx(i)+ajex(i)*ex(i)
            jy(i) = ajy(i)*jy(i)+ajey(i)*ey(i)
            jz(i) = ajz(i)*jz(i)+ajez(i)*ez(i)
        enddo
        !$omp end parallel do
            
end subroutine 