!-----------------------------------------------------------------------
!   完全導体
!-----------------------------------------------------------------------
subroutine PEC()   
   use constants 
   implicit none
   integer :: i,j
   real(kind=8) sp_rad,cy_rad
   write(30,*)"object location:"," ic:",ic," jc:",jc
   
  if(obj==1) then
    write(30,*)"object:cylinder"
    do j=0,ny
      do i=0,nx
        if(cy_rad(i,j)<=radius.and.cy_rad(i+1,j)<=radius) then
          aex(i,j)  =  0.0d0
          bexy(i,j) =  0.0d0
        endif
        if(cy_rad(i,j)<=radius.and.cy_rad(i,j+1)<=radius) then
          aey(i,j)  =  0.0d0
          beyx(i,j) =  0.0d0
        endif
        if(cy_rad(i,j)<=radius)then
          aez(i,j)    =  0.0d0
          bezx(i,j)   =  0.0d0
          bezy(i,j)   =  0.0d0
        endif
      end do
    end do
  else if(obj==2) then
    write(30,*)"object:proper prism"
    do j=jc-ly2,jc+ly2
      do i=ic-lx2,ic+lx2-1
        aex(i,j) = 0.0d0
        bexy(i,j) = 0.0d0
      enddo
    enddo
    do j=jc-ly2,jc+ly2-1
      do i=ic-lx2,ic+lx2
        aey(i,j) = 0.0d0
        beyx(i,j) = 0.0d0
      enddo
    enddo
    do j=jc-ly2,jc+ly2
      do i=ic-lx2,ic+lx2
        aez(i,j) = 0.0d0
        bezx(i,j) = 0.0d0
        bezy(i,j) = 0.0d0
      enddo
    enddo
  endif
end subroutine PEC