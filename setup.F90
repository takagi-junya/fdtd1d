subroutine setup()
    use HDF5
    use constants
    use hdfio
    implicit none
    real(kind=8) :: mux,muy,muz
    real(kind=8) :: a,epsx,epsy,epsz,sgmx,sgmy,sgmz,sgmm
    real(kind=8) :: sgex,sgey,sgez
    integer :: i,k
    real(kind=8) :: omg,ti,f

    namelist /space/nxx,dx,abc,lpml
    namelist /time/deltat,nstep
    namelist /output/interval,ostart,component,io
    namelist /wave/mode,amps,lambda,orgs,pt,pw
    namelist /plasma/pls,width,nu,wp,wc
        
    open(10,file="param.inp",action="read")
    read(10,nml=space)
    read(10,nml=time)
    read(10,nml=output)
    read(10,nml=wave)    
    read(10,nml=plasma)
    close(10)

    !吸収境界
    if(abc==0) then
        nx = nxx
    else
        nx = nxx+2*lpml(1)
    endif
    ic = int(0.5*nx)

    istart = int(nx/nprocs)*myrank+1 
    iend   = int(nx/nprocs)*(myrank+1)

    dimsf1d(1) = nx 
    dimsfi1d(1)= nx 
    dcount = iend-istart+1 
    chunk_dims1d(1) = dcount 

    allocate(ex(istart-1:iend+1),ey(istart-1:iend+1),ez(istart-1:iend+1))
    allocate(hx(istart-1:iend+1),hy(istart-1:iend+1),hz(istart-1:iend+1))
    allocate(jx(istart-1:iend+1),jy(istart-1:iend+1),jz(istart-1:iend+1))
    allocate(vx(istart-1:iend+1),vy(istart-1:iend+1),vz(istart-1:iend+1),nd(istart-1:iend+1))
    allocate(aex(istart-1:iend+1),aey(istart-1:iend+1),aez(istart-1:iend+1))
    allocate(bexy(istart-1:iend+1),beyx(istart-1:iend+1))
    allocate(bezx(istart-1:iend+1),bezy(istart-1:iend+1))
    allocate(amx(istart-1:iend+1),amy(istart-1:iend+1),amz(istart-1:iend+1))
    allocate(bmxy(istart-1:iend+1),bmyx(istart-1:iend+1))
    allocate(bmzx(istart-1:iend+1),bmzy(istart-1:iend+1))
    allocate(avx(istart-1:iend+1),avy(istart-1:iend+1),avz(istart-1:iend+1))
    allocate(epsd(istart-1:iend+1),sgmed(istart-1:iend+1))
    allocate(mud(istart-1:iend+1),sgmmd(istart-1:iend+1))
    allocate(ajx(istart-1:iend+1),ajy(istart-1:iend+1),ajz(istart-1:iend+1))
    allocate(ajex(istart-1:iend+1),ajey(istart-1:iend+1),ajez(istart-1:iend+1))
    allocate(omat(3,3),sa(3,3),sb(3,3),tc(3,3),sab(3,3),imat(3,3))

    filename(1:6)=(/"ex.h5","ey.h5","ez.h5","hx.h5","hy.h5","hz.h5"/)
    filename(7:10)=(/"wphi.h5","wz.h5","uphi.h5","uz.h5"/)
    filename(11:12)=(/"dphi.h5","dz.h5"/)
    filename(13:14)=(/"dphi.h5","dz.h5"/)
    filename(15) = "D.h5"
    filename(16:18)=(/"jx.h5","jy.h5","jz.h5"/)
    groupname(1:6)=(/"ex","ey","ez","hx","hy","hz"/)
    groupname(7:10)=(/"wphi","wz","uphi","uz"/)
    groupname(11:12)=(/"dphi","dz"/)
    groupname(13:14)=(/"phi","dz"/)
    groupname(15) = "D"
    groupname(16:18)=(/"jx","jy","jz"/)

    datasetname(1:6)=(/"wphi","wz","uphi","uz","dphi","dz"/)
    datasetname(7:9)=(/"dphi","dz","D"/)

    dt=deltat/(c*sqrt(1.0d0/(dx*dx)))
    freq = c/(lambda*dx)
    omega = 2*pai*freq
    k0 = omega/c
    jx = 0.0d0
    jy = 0.0d0
    jz = 0.0d0

    pc = orgs(1)*dx
    pw = pw*dx

    !背景媒質
    do i=istart,iend
        epsd(i)=epsbk
        mud(i)=mubk
        sgmed(i)=sigebk
        sgmmd(i)=sigmbk
    end do
    
    left = myrank-1
    if(myrank.eq.0) then
        left = mpi_proc_null 
    endif
    right = myrank+1
    if(myrank.eq.nprocs-1) then
        right = mpi_proc_null 
    endif

    !係数の計算
    call exchg1d(epsd,istart,iend,comm,left,right)
    call exchg1d(sgmed,istart,iend,comm,left,right)
    call exchg1d(mud,istart,iend,comm,left,right)
    call exchg1d(sgmmd,istart,iend,comm,left,right)

    !係数の計算
    do i=istart,iend
        epsx = 0.5d0*(epsd(i)+epsd(i))*eps0
        sgex = 0.5d0*(sgmed(i)+sgmed(i))
        a = 0.5d0*sgex*dt/epsx
        aex(i) = (1.0d0-a)/(1.0d0+a)

        epsy = 0.5d0*(epsd(i)+epsd(i-1))*eps0
        sgey = 0.5d0*(sgmed(i)+sgmed(i-1))
        a = 0.5d0*sgey*dt/epsy
        aey(i) = (1.0d0-a)/(1.0d0+a)
        beyx(i)= dt/epsy/(1.0d0+a)/dx
        
        epsz = 0.25*(epsd(i)+epsd(i-1)+epsd(i)+epsd(i-1))*eps0
        sgez = 0.25*(sgmed(i)+sgmed(i-1)+sgmed(i)+sgmed(i-1))
        a = 0.5d0*sgez*dt/epsz
        aez(i) = (1.0d0-a)/(1.0d0+a)
        bezx(i) = dt/epsz/(1.0d0+a)/dx 
        
        mux = 0.5d0*(mud(i)+mud(i-1))*mu0
        sgmx= 0.5d0*(sgmmd(i)+sgmmd(i-1))
        a = 0.5d0*sgmx*dt/mux
        amx(i) = (1.0d0-a)/(1.0d0+a)

        muy = 0.5d0*(mud(i)+mud(i))*mu0
        sgmy= 0.5d0*(sgmmd(i)+sgmmd(i))
        a = 0.5d0*sgmy*dt/muy
        amy(i) = (1.0d0-a)/(1.0d0+a)
        bmyx(i)= dt/muy/(1.0d0+a)/dx

        muz = mud(i)*mu0
        sgmz= sgmmd(i)
        a = 0.5d0*sgmz*dt/muz
        amz(i) = (1.0d0-a)/(1.0d0+a)
        bmzx(i) = dt/muz/(1.0d0+a)/dx
    enddo

    ipls = 0 
    iple = 0 

    if(istart.ge.ic-width.and.iend.le.ic+width) then
        ipls = istart
        iple = iend
    endif

    if(istart.le.ic-width.and.iend.le.ic+width) then
        ipls = ic-width 
        iple = iend 
    endif
    
    if(istart.le.ic-width.and.iend.ge.ic+width) then
        ipls = ic-width 
        iple = ic+width 
    endif

    if(istart.ge.ic-width.and.iend.ge.ic+width) then
        ipls = istart 
        iple = ic+width 
    endif

    if(pls.eq.1) then
        call ADE()
    elseif(pls.eq.2) then
        call JEC()
    elseif(pls.eq.3) then
        call EOM()
    endif

    if(myrank.eq.0) then
        open(30,file="sim.out")
        write(30,*)"nx:",nx
        write(30,'(a4,I4.4,a3,I4.4)')"ic:",ic
        write(30,*)"io:",io
        write(30,*)"frequency:",freq 
        write(30,*)"T:",1/freq
        write(30,*)"dt:",dt,"dx:",dx
    endif
    
    if(myrank.eq.nprocs-1) then 
        iend = nx-1
    endif
end subroutine setup