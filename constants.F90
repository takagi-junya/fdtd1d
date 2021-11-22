!--------------------------------------
!     シミュレーションで使用する定数、配列
!--------------------------------------
module constants
   use HDF5  
   use mpi            
    
   
   !全空間セル数
   integer :: nx

   !/space/
   integer :: nxx              !number of cells
   real(kind=8) :: dx          !grid length
   integer :: abc,lpml         !number of cells for absoption boundry

   !/time/
   real(kind=8) :: t
   real(kind=8) :: dt
   real(kind=8) :: deltat                    !CFL coefficent
   integer :: step,nstep                     !total time step

   !/output/
   integer :: interval                      !output interval 
   integer :: ostart                        !output start
   integer :: component(9)                  !EM-field component outputed
   integer :: io(2)

   !/object/
   integer :: ic
   
   !/wave/
   integer :: wshape
   real(kind=8) :: lambda,freq,omega      !wave length,frequency
   real(kind=8) :: amps(3)                !amplitude
   real(kind=8) :: orgs                      
   real(kind=8) :: pc,pt,pw,k0
   real(kind=8) :: alpha

   !/plasma
   integer :: pls
   real(kind=8) :: width !plasma slab length
   real(kind=8) :: wp,nu !plasma frequency,collision frequency
   real(kind=8) :: wc(3) !cyclotron frequency 
   real(kind=8),allocatable :: omat(:,:),sa(:,:),sb(:,:),tc(:,:),imat(:,:),sab(:,:)

   !PML
   integer,parameter :: order=4
   real(kind=8),parameter :: rmax=-120.0d0

   !HDF5
   integer :: h5count

   !MPI
   integer :: left,right,nprocs,myrank
   integer :: mpierr,comm,infom,hdferr
   integer :: istart,iend,ipls,iple,dcount
    
   !physical constants
   real(kind=8),parameter::eps0=8.854188d-12,mu0=1.256637d-6 !permitivity,permeability
   real(kind=8),parameter::qe=1.602176462d-19,mel = 9.10938188d-31 !charge,massive
   real(kind=8),parameter::c=2.9979246d8,z0=376.73031d0      !light speed,impedance
   real(kind=8),parameter::epsbk=1.0d0,mubk=1.0d0,sigebk=0.0d0,sigmbk=0.0d0
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