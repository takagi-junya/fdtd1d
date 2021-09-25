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

    namelist /space/nxx,nyy,dx,dy,abc,pbc,lpml
    namelist /time/deltat,nstep
    namelist /output/out,outs,comp,io,jo
    namelist /scatt/mode,lx,ly,gamma0,theta0,phi0,amp,lambda,tau0
    namelist /far/isx,isy,theta1,phi1
    namelist /object/obj,med,ic,jc,lx2,ly2,epsr,radius
    namelist /wave/kwave,amps,orgs,ang,pt,pw
    namelist /plasma/pls,prad,nu,wp,wc
        
    open(10,file="param.inp",action="read")
    read(10,nml=space)
    read(10,nml=time)
    read(10,nml=output)
    read(10,nml=scatt)
    read(10,nml=far)
    read(10,nml=object)
    read(10,nml=wave)    
    read(10,nml=plasma)
    close(10)

    !吸収境界
    if(abc==0) then
        nx = nxx
        ny = nyy
    else
        nx = nxx+2*lpml(1)
        ny = nyy+2*lpml(2)
    endif
    ic = int(0.5*nx)
    jc = int(0.5*ny)
    jo = jc
    write(30,*)"nx:",nx,"ny:",ny
    write(30,'(a3,I4.4,a3,I4.4)')"ic:",ic,"jc:"c

    allocate(ex(0:nx),ey(0:nx),ez(0:nx))
    allocate(hx(0:nx),hy(0:nx),hz(0:nx))
    allocate(jx(0:nx),jy(0:nx),jz(0:nx))
    allocate(vx(0:nx),vy(0:nx),vz(0:nx),nd(0:nx))
    allocate(aex(0:nx),aey(0:nx),aez(0:nx))
    allocate(bexy(0:nx),beyx(0:nx))
    allocate(bezx(0:nx),bezy(0:nx))
    allocate(amx(0:nx),amy(0:nx),amz(0:nx))
    allocate(bmxy(0:nx),bmyx(0:nx))
    allocate(bmzx(0:nx),bmzy(0:nx))
    allocate(avx(0:nx),avy(0:nx),avz(0:nx))
    allocate(epsd(-1:nx),sgmed(-1:nx))
    allocate(mud(-1:nx),sgmmd(-1:nx))
    allocate(ajx(0:nx),ajy(0:nx),ajz(0:nx))
    allocate(ajex(0:nx),ajey(0:nx),ajez(0:nx))
    

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

    dt=deltat/(c*sqrt(1.0d0/(dx*dx)+1.0d0/(dy*dy)))
    freq = c/(lambda*dx)
    omega = 2*pai*freq
    jx = 0.0d0
    jy = 0.0d0
    jz = 0.0d0

    if(kwave==0) then
        pc = orgs(1)*dx
        pw = pw*dx
    endif

    if(mode==1.or.mode==2.or.mode==3) then
        write(30,*)"frequency:",freq 
        write(30,*)"T:",1/freq
        write(30,*)"ka:",2*pai*(radius*dx)/(lambda*dx)
    endif

    write(30,*)"dt:",dt,"dx:",dx,"dy:",dy
    
    !背景媒質
    do i=-1,nx
        epsd(i)=epsbk
        mud(i)=mubk
        sgmed(i)=sigebk
        sgmmd(i)=sigmbk
    end do
    
    !係数の計算
    do i=0,nx
        epsx = 0.5d0*(epsd(i)+epsd(i))*eps0
        sgex = 0.5d0*(sgmed(i)+sgmed(i))
        a = 0.5d0*sgex*dt/epsx
        aex(i) = (1.0d0-a)/(1.0d0+a)
        bexy(i)= dt/epsx/(1.0d0+a)/dy

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
        bezy(i) = dt/epsz/(1.0d0+a)/dy
        
        mux = 0.5d0*(mud(i)+mud(i-1))*mu0
        sgmx= 0.5d0*(sgmmd(i)+sgmmd(i-1))
        a = 0.5d0*sgmx*dt/mux
        amx(i) = (1.0d0-a)/(1.0d0+a)
        bmxy(i)= dt/mux/(1.0d0+a)/dy

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
        bmzy(i) = dt/muz/(1.0d0+a)/dy
    enddo

    if(mode.eq.10) then
        write(30,*)"PLRC"
        allocate(phix(0:nx),phiy(0:nx),phiz(0:nx))
        allocate(aphix(0:nx),aphiy(0:nx),aphiz(0:nx))
        allocate(pex(0:nx),pey(0:nx),pez(0:nx))

        xhi0 = wp**2.0d0/nu*(dt-(1.0d0-exp(-nu*dt)/nu))
        xi0 = wp**2.0d0/nu*(dt**2.0d0/2.0d0+(dt+1.0d0/nu)*exp(-nu*dt)/nu-1.0d0/nu/nu)

        dxhi0 = -(wp/nu)**2.0d0*((exp(-nu*dt)-1.0d0)**2.0d0)
        dxi0 = (wp/nu)**2.0d0*(((dt+1.0d0/nu)*exp(-nu*dt)-1.0d0/nu)*(1.0d0-exp(-nu*dt)))
        do i=ic-prad,ic+prad+1
            aex(i) = (1.0d0-xi0)/(1.0d0+xhi0-xi0)
            bexy(i)= dt/eps0/(1.0d0+xhi0-xi0)/dy
            aphix(i) = 1.0d0/(1.0d0+xhi0-xi0)
        enddo

        do i=ic-prad,ic+prad
            aey(i) = (1.0d0-xi0)/(1.0d0+xhi0-xi0)
            beyx(i)= dt/eps0/(1.0d0+xhi0-xi0)/dx
            aphiy(i) = 1.0d0/(1.0d0+xhi0-xi0)
        
            aez(i) = (1.0d0-xi0)/(1.0d0+xhi0-xi0)
            bezx(i) = dt/eps0/(1.0d0+xhi0-xi0)/dx 
            bezy(i) = dt/eps0/(1.0d0+xhi0-xi0)/dy
            aphiz(i) = 1.0d0/(1.0d0+xhi0-xi0)
        enddo
        
    endif

    if(pls.eq.1) then
        call ADE()
        write(30,*)"ADE"
    elseif(pls.eq.2) then
        call JEC()
        write(30,*)"JEC"
    elseif(pls.eq.3) then
        call EOM()
        write(30,*)"EOM"
    endif

end subroutine setup