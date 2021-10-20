!-----------------------------------------------------------------------
!     PMLの共通変数 
!-----------------------------------------------------------------------
module pml2d !共通変数
   use constants
   implicit none
   type pml                                           
      integer::i0,i1                                        !PMLの範囲
      real(kind=8),pointer::expml(:),eypml(:),ezx(:),ezy(:) !PML電界
      real(kind=8),pointer::hxpml(:),hypml(:),hzx(:),hzy(:) !PML磁界
      real(kind=8),pointer::aexpml(:),aeypml(:)                 !諸係数
      real(kind=8),pointer::beypml(:),bexpml(:)
      real(kind=8),pointer::amxpml(:),amypml(:)
      real(kind=8),pointer::bmypml(:),bmxpml(:)
   end type pml
   type(pml)::pml_l,pml_r,pml_d,pml_u !左右，上下の４つPML領域の構造体
   real(kind=8)::copml
   parameter(copml=-1.5280063d-4)
   contains
  !-----------------------------------------------------------------------
  !     4PML領域初期設定のcall
  !-----------------------------------------------------------------------
   subroutine initpml()
      use constants
      implicit none
      if(myrank.eq.0) then
         write(30,*)"init_pml"
         write(30,*)"pmlx:",lpml(1)
         call init_pml(pml_l,0,lpml(1))                  !左側のPML
      endif
      if(myrank.eq.nprocs-1) then
         write(30,*)"pmlx2:",lpml(1)
         call init_pml(pml_r,nx-lpml(1),nx)              !右側のPML
      endif
   end subroutine initpml
  !-----------------------------------------------------------------------
  !     PMLの初期設定（配列の確保と係数の計算）
  !         p:PML構造体
  !         i0,i1:PMLの左辺，右辺の位置 
  !         j0,j1:PMLの下端，上端の位置 
  !-----------------------------------------------------------------------
   subroutine init_pml(p,i0,i1)
      use constants
      implicit none
      type(pml)::p                            !PML構造体
      integer,intent(in) :: i0,i1
      integer :: i
      real(kind=8)::a
      real(kind=8)::sigmxm,sigmxe
      real(kind=8)::sigmym,sigmye
      real(kind=8)::smax0x                      !x,y方向の導電率の最大値
      real(kind=8)::epspml,mupml                       !PMLの比誘電率，比透磁率
  
      !PMLの範囲を構造体に記憶
      p%i0=i0
      p%i1=i1
      !PML領域内の電磁界のスプリット成分の配列領域
      allocate(p%expml(i0:i1))   
      allocate(p%eypml(i0:i1))
      allocate(p%ezx(i0:i1))
      allocate(p%ezy(i0:i1))
      allocate(p%hxpml(i0:i1))
      allocate(p%hypml(i0:i1))
      allocate(p%hzx(i0:i1))
      allocate(p%hzy(i0:i1))
      !配列の初期化(平面波の場合はその初期値に設定）
      p%expml=0.0d0                       
      p%eypml=0.0d0
      p%ezx=0.0d0
      p%ezy=0.0d0
      p%hxpml=0.0d0
      p%hypml=0.0d0
      p%hzx=0.0d0
      p%hzy=0.0d0
      !係数の領域確保
      allocate(p%aeypml(i0:i1))
      allocate(p%aexpml(i0:i1))
      allocate(p%amypml(i0:i1))
      allocate(p%amxpml(i0:i1))
      allocate(p%beypml(i0:i1))
      allocate(p%bexpml(i0:i1))
      allocate(p%bmypml(i0:i1))
      allocate(p%bmxpml(i0:i1))
      !反射係数の要求精度〔dB〕から最大導電率を計算（真空）
      smax0x=copml*rmax*(order+1)/(lpml(1)*dx)
      !導電率分布，磁気伝導率分布
      mupml=mubk*mu0                                     !PMLの媒質定数
      epspml=epsbk*eps0
      do i=i0,i1        
         if(i<lpml(1)) then                              !左側のPML
            sigmxm=((lpml(1)-i-0.5d0)/lpml(1))**order*smax0x      
            sigmxe=(float(lpml(1)-i)/lpml(1))**order*smax0x
         else if(i>=nx-lpml(1)) then                     !右側のPML
            sigmxm=((i-nx+lpml(1)+0.5d0)/lpml(1))**order*smax0x
            sigmxe=(float(i-nx+lpml(1))/lpml(1))**order*smax0x
         else                                         
            sigmxm=0.0d0
            sigmxe=0.0d0
         end if
                                           
         sigmym=0.0d0
         sigmye=0.0d0
         !係数
         sigmxe=sigmxe*epsbk
         a=0.5d0*sigmxe*dt/epspml
         p%aexpml(i)=(1.0d0-a)/(1.0d0+a)
         p%bexpml(i)=dt/epspml/(1.0d0+a)/dx
         
         sigmye=sigmye*epsbk
         a=0.5d0*sigmye*dt/epspml
         p%aeypml(i)=(1.0d0-a)/(1.0d0+a)

         sigmxm=sigmxm*epsbk
         a=0.5d0*sigmxm*dt/epspml
         p%amxpml(i)=(1.0d0-a)/(1.0d0+a)
         p%bmxpml(i)=dt/mupml/(1.0d0+a)/dx

         sigmym=sigmym*epsbk
         a=0.5d0*sigmym*dt/epspml
         p%amypml(i)=(1.0d0-a)/(1.0d0+a)
      end do
   end subroutine init_pml
  !-----------------------------------------------------------------------
  !     4PML内の電界の計算
  !-----------------------------------------------------------------------
   subroutine epml()
      use constants
      implicit none
      if(myrank.eq.0) then
         call e_pml(pml_l)                        !左側のPML
      endif
      if(myrank.eq.nprocs-1) then
         call e_pml(pml_r)                        !右側のPML
      endif
   end subroutine epml
   !-----------------------------------------------------------------------
   !     電界の計算
   !-----------------------------------------------------------------------
   subroutine e_pml(p)
      
      use constants
      implicit none
      type(pml)::p
      integer::i,i0p,i1p

      i0p=p%i0                !PMLの領域：x方向：i0p to i1p
      i1p=p%i1
      !Ex
      do i=i0p,i1p-1 
         p%expml(i)=p%aeypml(i)*p%expml(i)
         ex(i)=p%expml(i)
      end do
      !Ey  
      do i=i0p+1,i1p-1
         p%eypml(i)=p%aexpml(i)*p%eypml(i)-p%bexpml(i)*(hz(i)-hz(i-1))
         ey(i)=p%eypml(i)
      end do
      !Ez
      do i=i0p+1,i1p-1 
         p%ezx(i)=p%aexpml(i)*p%ezx(i)+p%bexpml(i)*(hy(i)-hy(i-1))
         p%ezy(i)=p%aeypml(i)*p%ezy(i)
         ez(i)=p%ezx(i)+p%ezy(i)
      end do
   
   end subroutine e_pml
  !-----------------------------------------------------------------------
  !     4PML内の磁界の計算
  !-----------------------------------------------------------------------
   subroutine hpml()     
        
      use constants
      implicit none
      if(myrank.eq.0) then
         call h_pml(pml_l)                      !左側のPML
      endif
      if(myrank.eq.nprocs-1) then
         call h_pml(pml_r)                      !右側のPML
      endif
   end subroutine hpml
  !-----------------------------------------------------------------------
  !     磁界の計算
  !-----------------------------------------------------------------------
   subroutine h_pml(p)      
        
      use constants
      implicit none 
      type(pml)::p     
      integer::i,i0p,i1p
      
      i0p=p%i0                !PMLの領域：x方向：i0p to i1p
      i1p=p%i1
      !Hx     
      do i=i0p+1,i1p-1
         p%hxpml(i)=p%amypml(i)*p%hxpml(i)
         hx(i)=p%hxpml(i)
      end do
           
        
      !Hy   
      do i=i0p,i1p-1
         p%hypml(i)=p%amxpml(i)*p%hypml(i)+p%bmxpml(i)*(ez(i+1)-ez(i))
         hy(i)=p%hypml(i)
      end do
           
        
      !Hz
      do i=i0p,i1p-1     
         p%hzx(i)=p%amxpml(i)*p%hzx(i)-p%bmxpml(i)*(ey(i+1)-ey(i))
         p%hzy(i)=p%amypml(i)*p%hzy(i)
         hz(i)=p%hzx(i)+p%hzy(i)
      end do
   
   end subroutine h_pml
end module
  