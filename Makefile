FC		= ftn
CC		= cc

FFLAGS		= -ipo -O3 -no-prec-div -g -static -fopenmp -fp-model fast=2 -xMIC-AVX512 -fpp -Dssurf=3
CFLAGS		= -ipo -O3 -no-prec-div -g -static -fopenmp -fp-model fast=2 -xMIC-AVX512
LFLAGS		= -ipo -O3 -no-prec-div -g -static -fopenmp -fp-model fast=2 -xMIC-AVX512 -fpp -Dssurf=3

LDFLAGS		= 
HDF5DIR		= ./
FFTWDIR		= ./
MPIDIR = ./

LINKER		= $(FC)

MODS 	= constants.o  hdfio.o pwave.o RCPwave.o LCPwave.o pml1d.o output.o

FOBJS	= fdtd1d.o ADE.o JEC.o EOM.o velocity.o setup.o efield.o hfield.o current.o  currentEOM.o finalize.o comm1d.o

OBJS	= $(MODS) $(FOBJS) 

.SUFFIXES: .F90

fdtd:	$(OBJS)
		$(LINKER) $(OBJS) $(LFLAGS) $(LDFLAGS) -o $@

.F90.o:
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -I$(MPIDIR) -c $< -o $@
   
$(MODS):%.o:	%.F90 
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -c $< -o $@
$(KMODS):%.o:	%.F90
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -c $< -o $@
clean:;
		rm *.o; rm *.mod
