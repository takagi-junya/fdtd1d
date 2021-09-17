!-----------------------------------------------------------------------
!　　誘電佁E
!-----------------------------------------------------------------------
subroutine epsmu()   
   use constants
   implicit none
   real(kind=8) :: sp_rad,cy_rad
   integer :: i,j,k

   if(obj==1) then
     write(30,*)"object:cylinder"
       do j=0,ny
         do i=0,nx
           if(cy_rad(i,j)<=radius) then
             epsd(i,j)=epsr
             mud(i,j)=1.0d0
             sgmed(i,j)=0.0d0
             sgmmd(i,j)=0.0d0
           endif
         enddo
       enddo
   elseif(obj==2) then
     write(30,*)"object:quadrangular prism"
       do j=jc-ly2, jc+ly2-1
         do i=ic-lx2, ic+lx2-1
           epsd(i,j)=epsr
           mud(i,j)=1.0d0
           sgmed(i,j)=0.0d0
           sgmmd(i,j)=0.0d0
         end do
       end do
   endif
 
end subroutine epsmu
