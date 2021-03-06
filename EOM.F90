subroutine EOM()
    use constants 
    implicit none
    integer ::i,j
    integer ::ips,ipe
    ajj = dt/eps0 

    omat = dt/2.0d0*omat
    sa = adding_mat(imat(:,:),omat(:,:))
    sa = inversing_mat(sa(:,:))

    sb = adding_mat(imat(:,:),-omat(:,:))

    sab = matmul(sa,sb)
    tc = qe*dt/mel*sa
    
    do i=ipls,iple
        nd(i) = mel*eps0*(wp**2.0d0)/(qe**2.0d0)
        write(*,*)"nd:",nd(i)
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

    subroutine show_mat(x)
        use constants 
        implicit none 
        integer :: i,j 
        real(kind=8),intent(in) :: x(3,3)
        real(kind=8) :: y(3,3)
        do i = 1,3
            write(30,*),x(i,1),x(i,2),x(i,3)
        enddo
        write(30,*)""
    end subroutine


    function inversing_mat(x) result(ans)
        use constants
        implicit none
        integer :: n,lwork,ifail,info
        integer :: ipiv(3),work(192)
        integer :: i,j
        real(kind=8) :: ans(3,3),tmp(3,3)
        real(kind=8),intent(inout) :: x(3,3)
        n = 3
        lwork = n*64
        tmp = x
        call dgetrf(n,n,x,n,ipiv,info)
        call dgetri(n,x,n,ipiv,work,lwork,info)
        tmp = matmul(tmp,x)
        ans = x
    end function
        
end subroutine
            