!------------------------------------------------------
!     constants and arrays 
!     auther takagi junya
!------------------------------------------------------
module constants
   use HDF5              
    
   !/space/
   integer :: nxx,nyy              !number of cells
   real(kind=8) :: dx,dy           !grid length
   integer :: abc                  !absorption boundry flag
   integer :: pbc(3)               !periodic boundry flag
   integer :: lpml(3)              !number of absorption boundry cells

   !/time/
   real(kind=8) :: t
   real(kind=8) :: dt
   real(kind=8) :: deltat          !CFL coefficient
   integer :: step,nstep           !number of time steps

   !/output/
   integer :: out                      !output interval
   integer :: outs                     !first output time
   integer :: comp(9)                  !output conponent
   integer :: io(2),jo(2)              !two output points
               
   !/scat/
   integer :: mode                     
   integer ::lx,ly                     !scatterd field cells
   real(kind=8) :: gamma0,theta0,phi0,amp !pwave angles,amplitude
   real(kind=8) :: lambda,tau0,freq,omega !wave length,time constant frequancy,angle frequancy
   real(kind=8),allocatable :: omat(:,:),omata(:,:),omatb(:,:),omatc(:,:),imat(:,:) !matrix for EOM

   !/far/
   integer :: isx,isy
   real(kind=8) :: theta1,phi1
   
   !/object/
   integer :: obj                      !1:cylinder 2:prism
   integer :: med                      !1:dielectric 2:PEC 
   integer :: ic,jc                    !center of simulation filed
   integer :: lx2,ly2                  !length of prism edge
   real(kind=8) :: epsr,radius         !relative permittivity cylinder's radius
   
   !/feed/
   integer :: ip,jp,kp
   real(kind=8) :: duration,t0
   
   !/wave/
   integer :: kwave
   real(kind=8) :: tau
   real(kind=8) :: amps(6)       !pwave angle
   real(kind=8) :: orgs(3)    
   real(kind=8) :: ang(3)        !pwave angle      
   real(kind=8) :: pc,pt,pw
   real(kind=8) :: alpha,k0

   !/plasma
   integer :: pls
   real(kind=8) :: prad,nu
   real(kind=8) :: wp
   real(kind=8) :: wc(3)

   integer ::ifed,jfed,kfed
   !PML
   integer,parameter :: order=4
   real(kind=8),parameter :: rmax=-120.0d0

   integer :: nx,ny
   
   real(kind=8),parameter::epsbk=1.0d0,mubk=1.0d0,sigebk=0.0d0,sigmbk=0.0d0

   !HDF5
   integer :: h5count

    
   real(kind=8),parameter::eps0=8.854188d-12,mu0=1.256637d-6 !permitivity and permeability of vacuum
   real(kind=8),parameter::qe=1.602176462d-19,mel = 9.10938188d-31 !quantity of charge,mass
   real(kind=8),parameter::c=2.9979246d8,z0=376.73031d0      !speed of light ,impedance
   real(kind=8),parameter::radi0=1.74532925d-2              
   real(kind=8),parameter::pai=3.14159265d0                     

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