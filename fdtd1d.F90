!-----------------------------
!   1次元FDTD
!-----------------------------
program fdtd1d
    use constants
    use HDF5
    use mpi
    use hdfio
    use pml2d
    use pwave
    use RCPwave
    use LCPwave
    use omp_lib
    use output
    implicit none
    real(kind=8) :: time0,time1
    integer :: nthread

    comm = mpi_comm_world 
    infom = mpi_info_null 
    call mpi_init(mpierr)
    call mpi_comm_size(comm,nprocs,mpierr)
    call mpi_comm_rank(comm,myrank,mpierr)
    
    !$omp parallel
    nthread = omp_get_num_threads()
    !$omp end parallel

    call hdfinit()
    call setup()
    call output_init()

    if((wshape.eq.0).or.(wshape.eq.1)) then
        call init_pwave()
    elseif(wshape.eq.2) then
        call RCP()
    elseif(wshape.eq.3) then
        call LCP()
    endif
    call initpml()

    t=dt

    call mpi_barrier(comm,mpierr)
    time0 = mpi_wtime()
    do step=1,nstep+1
        if(myrank.eq.0) then
            write(*,'(a11,I7.4)')"Time step:",step
        endif
        call efield()
        call epml()
        t=t+0.5d0*dt
        call hfield()
        if(pls.eq.1) then
            call current()
        endif
        if(pls.eq.2) then
            call current()
        endif
        if(pls.eq.3) then
            call velocity()
            call currentEOM()
        endif
        call hpml()
        t=t+0.5d0*dt
        call out_emf(step-1)
    enddo
    call mpi_barrier(comm,mpierr)
    time1 = mpi_wtime()

    if(myrank.eq.0) then
        open(50,file="memo.txt",position="append")
        write(50,*)"np:",nprocs,"nt:",nthread," time:",time1-time0
        close(50)
    endif

    call hdffinalize()
    call finalize()
    call output_fin()
    call mpi_finalize(mpierr)
end program fdtd1d
