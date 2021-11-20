
!-----------------------------------------------------------------------
!       計算結果の出力
!-----------------------------------------------------------------------
module output
   integer :: is,ie
   contains
   subroutine output_init()
      use HDF5
      use mpi 
      use constants
      use hdfio
      implicit none

      is=lpml
      ie=nx-lpml
      dims(1)=ie-is+1

      if(myrank.eq.0) then
         write(30,*)"io1:",io(1)
         write(30,*)"io2:",io(2)
      endif

      if(pls.eq.0) then
         io(1) = ic

         if(io(1).ge.istart.and.io(1).le.iend) then
            open(101,file="exin.txt")
            open(102,file="eyin.txt")
            open(103,file="ezin.txt")
         endif
      
      else 
         if(io(1).ge.istart.and.io(1).lt.iend) then
            write(*,*)"io1:",io(1)
            open(104,file="ex_ref.txt")
            open(105,file="ey_ref.txt")
            open(106,file="ez_ref.txt")
         endif
         call mpi_barrier(comm,mpierr)
         if(io(2).ge.istart.and.io(2).lt.iend) then
            write(*,*)"io2:",io(2)
            open(107,file="ex_trs.txt")
            open(108,file="ey_trs.txt")
            open(109,file="ez_trs.txt")
         endif
      endif

      if(component(1).eq.1) then
         call h5popen(filename(1),groupname(1),file_id(1),group_id(1),plist_id(1),0,comm,infom,hdferr)
      endif
      if(component(2).eq.1) then
         call h5popen(filename(2),groupname(2),file_id(2),group_id(2),plist_id(2),0,comm,infom,hdferr)
      endif
      if(component(3).eq.1) then
         call h5popen(filename(3),groupname(3),file_id(3),group_id(3),plist_id(3),0,comm,infom,hdferr)
      endif
      if(component(4).eq.1) then
         call h5popen(filename(4),groupname(4),file_id(4),group_id(4),plist_id(4),0,comm,infom,hdferr)
      endif
      if(component(5).eq.1) then
         call h5popen(filename(5),groupname(5),file_id(5),group_id(5),plist_id(5),0,comm,infom,hdferr)
      endif
      if(component(6).eq.1) then
         call h5popen(filename(6),groupname(6),file_id(6),group_id(6),plist_id(6),0,comm,infom,hdferr)
      endif
      if(component(7).eq.1) then
         call h5popen(filename(16),groupname(16),file_id(7),group_id(7),plist_id(7),0,comm,infom,hdferr)
      endif
      if(component(8).eq.1) then
         call h5popen(filename(17),groupname(17),file_id(8),group_id(8),plist_id(8),0,comm,infom,hdferr)
      endif
      if(component(9).eq.1) then
         call h5popen(filename(18),groupname(18),file_id(9),group_id(9),plist_id(9),0,comm,infom,hdferr)
      endif

   end subroutine 

   subroutine out_emf(n)
      use HDF5
      use mpi 
      use constants
      use hdfio
      implicit none
   
      integer :: i,j,k,l
      integer,intent(in) :: n

      if(pls.eq.0) then
         io(1) = ic

         if(io(1).ge.istart.and.io(1).le.iend) then
            write(101,*)ex(io(1))
            write(102,*)ey(io(1))
            write(103,*)ez(io(1))
         endif
      
      else 
         if(io(1).ge.istart.and.io(1).le.iend) then
            write(104,*)ex(io(1))
            write(105,*)ey(io(1))
            write(106,*)ez(io(1))
         endif

         if(io(2).ge.istart.and.io(2).le.iend) then
            write(107,*)ex(io(2))
            write(108,*)ey(io(2))
            write(109,*)ez(io(2))
         endif
      endif
      if(mod(int(n-ostart),interval).eq.0) then  
         write(tag,'(I4.4)') h5count
         if(myrank.eq.0) then
            write(*,'(/a8,I6.5/)')'h5 out:',n
         endif
         if(component(1).eq.1) then
            call wrt1p(group_id(1),plist_id(1),tag,dimsf1d,chunk_dims1d,myrank,ex(istart:iend))
         endif
         if(component(2).eq.1) then
            call wrt1p(group_id(2),plist_id(2),tag,dimsf1d,chunk_dims1d,myrank,ey(istart:iend))
         endif
         if(component(3).eq.1) then
            call wrt1p(group_id(3),plist_id(3),tag,dimsf1d,chunk_dims1d,myrank,ez(istart:iend))
         endif
         if(component(4).eq.1) then
            call wrt1p(group_id(4),plist_id(4),tag,dimsf1d,chunk_dims1d,myrank,hx(istart:iend))
         endif
         if(component(5).eq.1) then
            call wrt1p(group_id(5),plist_id(5),tag,dimsf1d,chunk_dims1d,myrank,hy(istart:iend))
         endif
         if(component(6).eq.1) then
            call wrt1p(group_id(6),plist_id(6),tag,dimsf1d,chunk_dims1d,myrank,hz(istart:iend))
         endif
         if(component(7).eq.1) then
            call wrt1p(group_id(7),plist_id(7),tag,dimsf1d,chunk_dims1d,myrank,jx(istart:iend))
         endif
         if(component(8).eq.1) then
            call wrt1p(group_id(8),plist_id(8),tag,dimsf1d,chunk_dims1d,myrank,jy(istart:iend))
         endif
         if(component(9).eq.1) then
            call wrt1p(group_id(9),plist_id(9),tag,dimsf1d,chunk_dims1d,myrank,jz(istart:iend))
         endif
	      h5count = h5count + 1
      endif 
   end subroutine

   subroutine output_fin()
      use constants
      use HDF5 
      use mpi 
      use hdfio 
      
      if(pls.eq.0) then
         io(1) = ic

         if(io(1).ge.istart.and.io(1).le.iend) then
            close(101)
            close(102)
            close(103)
         endif
      
      else 
         if(io(1).ge.istart.and.io(1).lt.iend) then
            close(104)
            close(105)
            close(106)
         endif

         if(io(2).ge.istart.and.io(2).lt.iend) then
            close(107)
            close(108)
            close(109)
         endif
      endif

      if(component(1).eq.1) then
         call h5pclose(file_id(1),group_id(1),plist_id(1),istat1(1))
      endif
      if(component(2).eq.1) then
         call h5pclose(file_id(2),group_id(2),plist_id(2),istat1(2))
      endif
      if(component(3).eq.1) then
         call h5pclose(file_id(3),group_id(3),plist_id(3),istat1(3))
      endif
      if(component(4).eq.1) then
         call h5pclose(file_id(4),group_id(4),plist_id(4),istat1(4))
      endif
      if(component(5).eq.1) then
         call h5pclose(file_id(5),group_id(5),plist_id(5),istat1(5))
      endif
      if(component(6).eq.1) then
         call h5pclose(file_id(6),group_id(6),plist_id(6),istat1(6))
      endif
      if(component(7).eq.1) then
         call h5pclose(file_id(7),group_id(7),plist_id(7),istat1(7))
      endif
      if(component(8).eq.1) then
         call h5pclose(file_id(8),group_id(8),plist_id(8),istat1(8))
      endif
      if(component(9).eq.1) then
         call h5pclose(file_id(9),group_id(9),plist_id(9),istat1(9))
      endif
   end subroutine
end module
