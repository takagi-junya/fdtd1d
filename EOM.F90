!------------------------------------------------------
!   Equation of motion for drude plasma model
!------------------------------------------------------
subroutine EOM()
    use constants 
    implicit none
    integer ::i,j
    real(kind=8) add_mat(3,3),inverse_mat(3,3)
    allocate(omat(3,3),omata(3,3),omatb(3,3),omatc(3,3),imat(3,3))
    ajj = dt/eps0 

    imat(1,1) = 1
    imat(2,2) = 1
    imat(3,3) = 1

    omat(1,1) = nu
    omat(2,2) = nu 
    omat(3,3) = nu 
    omat(1,2) = wc(3)
    omat(1,3) = -wc(2)
    omat(2,1) = -wc(3)
    omat(2,3) = wc(1)
    omat(3,1) = wc(2)
    omat(3,2) = -wc(1)
    omat = dt/2.0d0*omat

    omata = adding_mat(imat(:,:),omat(:,:))
    omata = inverseing_mat(omata(:,:))
    omatb = adding_mat(imat(:,:),-omat(:,:))
    omatc = qe*dt/mel*omata

    omata = matmul(omata,omatb)
    omatb = omatc
    do i=ic-prad,ic+prad-1
        avx(i) = 1.0d0
        ajex(i) = -(2.0d0*qe/mel*dt)/(2.0d0+nu*dt)
    enddo

    do i=ic-prad,ic+prad
        avy(i) = 1.0d0
        ajey(i) = -(2.0d0*qe/mel*dt)/(2.0d0+nu*dt)
    enddo

    do i=ic-prad,ic+prad
        avz(i) = 1.0d0
        ajez(i) = -(2.0d0*qe/mel*dt)/(2.0d0+nu*dt)
        nd(i) = mel*eps0*(wp**2.0d0)/(qe**2.0d0)
    enddo

    contains
        
    function adding_mat(x,y) result(ans)
        use constants
        implicit none 
        integer :: i,j
        real(kind=8),intent(in) :: x(3,3),y(3,3)
        real(kind=8) :: ans(3,3)
        
        do i=1,3
            do j=1,3
                ans(i,j) = x(i,j) + y(i,j)
            enddo
        enddo
        
    end function 


    function inverseing_mat(x) result(ans)
        use constants
        implicit none
        integer :: n,lwork,ifail,info
        integer :: ipiv(3),work(192)
        integer :: i,j
        real(kind=8) :: ans(3,3)
        real(kind=8),intent(inout) :: x(3,3)
        n = 3
        lwork = n*64
        call dgetrf(n,n,x,n,ipiv,info)
        call dgetri(n,x,n,ipiv,work,lwork,info)
        ans = x
    end function
        
end subroutine
            