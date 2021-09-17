real(kind=8) function pscat(x)
    use constants
    implicit none
    real(kind=8),intent(in) :: x
    real(kind=8) :: xx,a
    xx = x-c*t-pc
    a = xx/pw
    pscat = amp*exp(-a*a)
end function pscat

real(kind=8) function pscat2(x)
    use constants
    implicit none
    real(kind=8),intent(in) :: x
    real(kind=8) :: xx,a
    xx = x-c*t-pc
    a = xx/pw
    pscat2 = amp*exp(-a*a)*cos(2.0d0*pai*freq/c*xx)
end function pscat2


real(kind=8) function sp_rad(x,y)
   use constants
   implicit none
   integer,intent(in) :: x,y
   real(kind=8) :: xx,yy
   xx = x
   yy = y
   sp_rad = sqrt((xx-ic)**2+(yy-jc)**2)
end function sp_rad

real(kind=8) function cy_rad(x,y)
   use constants
   implicit none
   integer,intent(in) :: x,y
   real(kind=8) :: xx,yy
   xx = x
   yy = y
   cy_rad = sqrt((xx-ic)**2+(yy-jc)**2)
end function cy_rad

function add_mat2(x,y)
    use constants
    implicit none 
    integer :: i,j
    real(kind=8),intent(in) :: x(3,3),y(3,3)
    real(kind=8) :: add_mat2(3,3),result(3,3)
    do i=1,3
        do j=1,3
            result(i,j) = x(i,j) + y(i,j)
        enddo
    enddo
    add_mat2 = result
end function 
