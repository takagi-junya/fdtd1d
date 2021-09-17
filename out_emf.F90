
!-----------------------------------------------------------------------
!       計算結果の出劁E
!-----------------------------------------------------------------------
subroutine out_emf(n)
   use HDF5
   use constants
   use hdfio
   implicit none
  
   integer :: i,j,k,l
   integer,intent(in) :: n
   integer :: is,ie

   is=lpml(1)
   ie=nx-lpml(1)
   dims(1)=ie-is+1
   
   if(n==0) then
     write(30,*)"io:",io
   endif
   open(100,file="ex.txt")
   write(100,*)ex(io)
   if(obj==0) then
      open(101,file="eyin.txt")
      write(101,*)ey(io)
   else
      open(103,file="eyinc.txt")
      write(103,*)ey(io)
   endif
   open(102,file="ez.txt")
   write(102,*)ez(io)
   
   if((n-outs)==0) then
      if(comp(1)==1) then
         call hdfopen(filename(1),groupname(1),file_id(1),group_id(1),0)
         call wrt1d(file_id(1),group_id(1),"0000",dims,ex(is:ie),istat1(1),istat2(1))
         call hdfclose(file_id(1),group_id(1),error(1))
      endif
      if(comp(2)==1) then
         call hdfopen(filename(2),groupname(2),file_id(2),group_id(2),0)
         call wrt1d(file_id(2),group_id(2),"0000",dims,ey(is:ie),istat1(2),istat2(2))
         call hdfclose(file_id(2),group_id(2),error(2))
      endif
      if(comp(3)==1) then
         call hdfopen(filename(3),groupname(3),file_id(3),group_id(3),0)
         call wrt1d(file_id(3),group_id(3),"0000",dims,ez(is:ie),istat1(3),istat2(3))
         call hdfclose(file_id(3),group_id(3),error(3))
      endif
      if(comp(4)==1) then
         call hdfopen(filename(4),groupname(4),file_id(4),group_id(4),0)
         call wrt1d(file_id(4),group_id(4),"0000",dims,hx(is:ie),istat1(4),istat2(4))
         call hdfclose(file_id(4),group_id(4),error(4))
      endif
      if(comp(5)==1) then
         call hdfopen(filename(5),groupname(5),file_id(5),group_id(5),0)
         call wrt1d(file_id(5),group_id(5),"0000",dims,hy(is:ie),istat1(5),istat2(5))
         call hdfclose(file_id(5),group_id(5),error(5))
      endif
      if(comp(6)==1) then
         call hdfopen(filename(6),groupname(6),file_id(6),group_id(6),0)
         call wrt1d(file_id(6),group_id(6),"0000",dims,hz(is:ie),istat1(6),istat2(6))
         call hdfclose(file_id(6),group_id(6),error(6))
      endif
      if(comp(7)==1) then
         call hdfopen(filename(16),groupname(16),file_id(16),group_id(16),0)
         call wrt1d(file_id(16),group_id(16),"0000",dims,jx(is:ie),istat1(16),istat2(16))
         call hdfclose(file_id(16),group_id(16),error(16))
      endif
      if(comp(8)==1) then
         call hdfopen(filename(17),groupname(17),file_id(17),group_id(17),0)
         call wrt1d(file_id(17),group_id(17),"0000",dims,jy(is:ie),istat1(17),istat2(17))
         call hdfclose(file_id(17),group_id(17),error(17))
      endif
      if(comp(9)==1) then
         call hdfopen(filename(18),groupname(18),file_id(18),group_id(18),0)
         call wrt1d(file_id(18),group_id(18),"0000",dims,jz(is:ie),istat1(18),istat2(18))
         call hdfclose(file_id(18),group_id(18),error(18))
      endif
   elseif(mod(int(n-outs),out)==0.and.n-outs>0.0d0) then
      h5count = h5count + 1  
      write(tag,'(I4.4)') h5count
      write(*,'(/,a15,a4,/)')'h5 out:',tag
      if(comp(1)==1) then
         call hdfopen(filename(1),groupname(1),file_id(1),group_id(1),1)
         call wrt1d(file_id(1),group_id(1),tag,dims,ex(is:ie),istat1(1),istat2(1))
         call hdfclose(file_id(1),group_id(1),error(1))
      endif
      if(comp(2)==1) then
         call hdfopen(filename(2),groupname(2),file_id(2),group_id(2),1)
         call wrt1d(file_id(2),group_id(2),tag,dims,ey(is:ie),istat1(2),istat2(2))
         call hdfclose(file_id(2),group_id(2),error(2))
      endif
      if(comp(3)==1) then
         call hdfopen(filename(3),groupname(3),file_id(3),group_id(3),1)
         call wrt1d(file_id(3),group_id(3),tag,dims,ez(is:ie),istat1(3),istat2(3))
         call hdfclose(file_id(3),group_id(3),error(3))
      endif
      if(comp(4)==1) then
         call hdfopen(filename(4),groupname(4),file_id(4),group_id(4),1)
         call wrt1d(file_id(4),group_id(4),tag,dims,hx(is:ie),istat1(4),istat2(4))
         call hdfclose(file_id(4),group_id(4),error(4))
      endif
      if(comp(5)==1) then
         call hdfopen(filename(5),groupname(5),file_id(5),group_id(5),1)
         call wrt1d(file_id(5),group_id(5),tag,dims,hy(is:ie),istat1(5),istat2(5))
         call hdfclose(file_id(5),group_id(5),error(5))
      endif
      if(comp(6)==1) then
         call hdfopen(filename(6),groupname(6),file_id(6),group_id(6),1)
         call wrt1d(file_id(6),group_id(6),tag,dims,hz(is:ie),istat1(6),istat2(6))
         call hdfclose(file_id(6),group_id(6),error(6))
      endif
      if(comp(7)==1) then
         call hdfopen(filename(16),groupname(16),file_id(16),group_id(16),1)
         call wrt1d(file_id(16),group_id(16),tag,dims,jx(is:ie),istat1(16),istat2(16))
         call hdfclose(file_id(16),group_id(16),error(16))
      endif
      if(comp(8)==1) then
         call hdfopen(filename(17),groupname(17),file_id(17),group_id(17),1)
         call wrt1d(file_id(17),group_id(17),tag,dims,jy(is:ie),istat1(17),istat2(17))
         call hdfclose(file_id(17),group_id(17),error(17))
      endif
      if(comp(8)==1) then
         call hdfopen(filename(18),groupname(18),file_id(18),group_id(18),1)
         call wrt1d(file_id(18),group_id(18),tag,dims,jz(is:ie),istat1(18),istat2(18))
         call hdfclose(file_id(18),group_id(18),error(18))
      endif
   endif 
end subroutine out_emf