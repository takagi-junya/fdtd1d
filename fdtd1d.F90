program fdtd2d
    use constants
    use HDF5
    use hdfio
    use pml2d
    use pwave
    use RCP 
    use LCP
    use omp_lib
    implicit none
    open(30,file="sim.out")
    call hdfinit()
    call setup()
    if(mode.le.1) then
        write(30,*)"set Plane wave"
        call init_pwave()
        call initpml()
    else if(mode.eq.2) then
        write(30,*)"set RCP wave"
        call setRCP()
        call initpml()
    else if(mode.eq.3) then
        write(30,*)"set LCP wave"
        call setLCP()
        call initpml()
    endif
    t=dt
    if(mode.eq.0) then
        write(30,*),"All total field mode"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call epml()
            t=t+0.5d0*dt
            call hfield()
            call current()
            call hpml()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
    elseif(mode.ge.1) then
        write(30,*),"All total field mode"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call epml()
            t=t+0.5d0*dt
            call velocity()
            call currentEOM()
            call hfield()
            call hpml()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
    endif
    call hdffinalize()
    close(30)
end program fdtd2d