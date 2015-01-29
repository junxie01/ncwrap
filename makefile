## makefile for netcdf example programs

# initialize suffix
.SUFFIXES:
.SUFFIXES:.F90 .o .F90

FC     = gfortran # Compiler
NCLIB  = /usr/local/lib # netcdf library path
NCINC  = /usr/local/include # netcdf module file (netcdf.mod)
FFLAGS =  #
NCFLAGS = -D_NETCDF -lnetcdf -lnetcdff #( -lnetcdff is required for netcdf ver 4)

#FC     = ifort # Compiler
#NCLIB  = ${HOME}/local/lib # netcdf library path
#NCINC  = ${HOME}/local/include # netcdf module file (netcdf.mod)
#FFLAGS =  #
#NCFLAGS = -D_NETCDF -lnetcdf #-lnetcdff #( -lnetcdff is required for netcdf ver 4)

all: test__netcdf_write.x test__netcdf_read.x

m_ncwrap.o: m_ncwrap.F90
	$(FC) -I$(NCINC) -L$(NCLIB) $(FFLAGS) $(NCFLAGS) -c $< 
test__netcdf_write.x: test__netcdf_write.F90 m_ncwrap.o
	$(FC) -I$(NCINC) -L$(NCLIB) $(FFLAGS) $(NCFLAGS) $^ -o $@
test__netcdf_read.x: test__netcdf_read.F90 m_ncwrap.o
	$(FC) -I$(NCINC) -L$(NCLIB) $(FFLAGS) $(NCFLAGS) $^ -o $@

clean:
	/bin/rm test__netcdf_write.x test__netcdf_read.x *.o *.mod
