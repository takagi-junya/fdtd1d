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
    if(mode.eq.0.or.mode.eq.1) then
        call init_pwave()
        call initpml()
    else if(mode.eq.2) then
        call RCP()
        call initpml()
    else if(mode.eq.3) then
        call LCP()
        call initpml()
    endif
    t=dt
    if(mode.eq.0) then
        write(30,*),"All total field mode"
        do step=1,nstep+1
            write(*,'(a10,I6.6)')"Time step:",step-1
            call efield()
            call epml()
            t=t+0.5d0*dt
            call hfield()
            call current()
            call hpml()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
    elseif((mode.eq.1).or.(mode.eq.2).or.(mode.eq.3)) then
        write(30,*),"All total field mode"
        call mpi_barrier(comm,mpierr)
        time0 = mpi_wtime()
        do step=1,nstep+1
            !write(*,'(a10,I6.6)')"Time step:",step-1
            call efield()
            call epml()
            t=t+0.5d0*dt
            call hfield()
            call velocity
            call currentEOM()
            call hpml()
            t=t+0.5d0*dt
            call out_emf(step-1)
            !call mpi_barrier(comm,mpierr)
        enddo
        call mpi_barrier(comm,mpierr)
        time1 = mpi_wtime()
    endif
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
