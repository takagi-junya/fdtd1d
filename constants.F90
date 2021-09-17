module constants
   use HDF5              
    
   !/space/
   integer :: nxx,nyy              !計算領域セル数
   real(kind=8) :: dx,dy                   !グリッド幅
   integer :: abc
   integer :: pbc(3),lpml(3) !吸収境界セル数

   !/time/
   real(kind=8) :: t
   real(kind=8) :: dt
   real(kind=8) :: deltat                      !CFL係数
   integer :: step,nstep                     !総step数

   !/output/
   integer :: out                      !出力間隔
   integer :: outs                     
   integer :: comp(9)                  !出力電磁場
   integer :: io,jo  
                  !出力ポイント
   !/scat/
   integer :: mode                     !0:FDTD 1:TF-SF 2:遠方界計算 
   integer ::lx,ly                     !散乱界領域セル数
   real(kind=8) :: gamma0,theta0,phi0,amp
   real(kind=8) :: lambda,tau0,freq,omega
   real(kind=8) :: bf(3)
   real(kind=8),allocatable :: omat(:,:),omata(:,:),omatb(:,:),omatc(:,:),imat(:,:)

   !/far/
   integer :: isx,isy
   real(kind=8) :: theta1,phi1
   
   !/object/
   integer :: obj                      !1:円柱 2:プリズム
   integer :: med                      !1:誘電体2:完全導体  
   integer :: ic,jc                    !散乱体の中心  
   integer :: lx2,ly2                  !プリズムの辺の長さ   
   real(kind=8) :: epsr,radius                 !比誘電率,半径   
   
   !/feed/
   integer :: ip,jp,kp
   real(kind=8) :: duration,t0
   
   !/wave/
   integer :: kwave
   real(kind=8) :: tau
   real(kind=8) :: amps(6)
   real(kind=8) :: orgs(3)
   real(kind=8) :: ang(3)              
   real(kind=8) :: pc,pt,pw
   real(kind=8) :: alpha

   !/plasma
   integer :: pls
   real(kind=8) :: prad,nu
   real(kind=8) :: wp
   real(kind=8) :: wc(3)

   integer ::ifed,jfed,kfed
   !PML
   integer,parameter :: order=4
   real(kind=8),parameter :: rmax=-120.0d0

   !全空間セル数
   integer :: nx,ny
   
   real(kind=8),parameter::epsbk=1.0d0,mubk=1.0d0,sigebk=0.0d0,sigmbk=0.0d0

   !HDF5
   integer :: h5count

    
   real(kind=8),parameter::eps0=8.854188d-12,mu0=1.256637d-6 !真空の誘電玁E��E��磁率
   real(kind=8),parameter::qe=1.602176462d-19,mel = 9.10938188d-31
   real(kind=8),parameter::c=2.9979246d8,z0=376.73031d0      !光速、インピ�Eダンス
   real(kind=8),parameter::radi0=1.74532925d-2               !角度のラジアン変換
   real(kind=8),parameter::pai=3.14159265d0                  !冁E��玁E   

   real(kind=8),allocatable :: ex(:),ey(:),ez(:)
   real(kind=8),allocatable :: pex(:),pey(:),pez(:)
   real(kind=8),allocatable :: hx(:),hy(:),hz(:)
   real(kind=8),allocatable :: vx(:),vy(:),vz(:)
   real(kind=8),allocatable :: jx(:),jy(:),jz(:),nd(:)
    
   real(kind=8),allocatable :: aex(:),aey(:),aez(:)
   real(kind=8),allocatable :: bexy(:),bexz(:)
   real(kind=8),allocatable :: beyx(:),beyz(:)
   real(kind=8),allocatable :: bezx(:),bezy(:)

   real(kind=8),allocatable :: amx(:),amy(:),amz(:)
   real(kind=8),allocatable :: bmxy(:),bmxz(:)
   real(kind=8),allocatable :: bmyx(:),bmyz(:)
   real(kind=8),allocatable :: bmzx(:),bmzy(:)

   real(kind=8),allocatable :: avx(:),avy(:),avz(:)
   real(kind=8),allocatable :: epsd(:),sgmed(:),mud(:),sgmmd(:)

   real(kind=8) :: ajj
   real(kind=8),allocatable :: ajx(:),ajy(:),ajz(:)
   real(kind=8),allocatable :: ajex(:),ajey(:),ajez(:)

   real(kind=8),allocatable :: phix(:),phiy(:),phiz(:)
   real(kind=8),allocatable :: aphix(:),aphiy(:),aphiz(:) 
   real(kind=8) :: xi0,xhi0,dxi0,dxhi0

    
end module constants