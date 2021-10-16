subroutine exchg1d(data,s,e,comm1d,leftp,rightp)
    use constants
    use mpi
    implicit none
    integer,intent(in) :: s,e,comm1d,leftp,rightp
    real(kind=8),intent(in) :: data(s-1:e+1)
    call MPI_SENDRECV(data(e),1,MPI_DOUBLE_PRECISION,rightp,0,&
                     &data(s-1),1,MPI_DOUBLE_PRECISION,leftp,0,comm1d,mpi_status_ignore,mpierr)
    call MPI_SENDRECV(data(s),1,MPI_DOUBLE_PRECISION,leftp,1,&
                    &data(e+1),1,MPI_DOUBLE_PRECISION,rightp,1,comm1d,mpi_status_ignore,mpierr)
end subroutine 