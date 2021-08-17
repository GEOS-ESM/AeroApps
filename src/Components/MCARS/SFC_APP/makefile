BASEDIR = /discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v3.2.0_buildSLES11SP1/x86_64-unknown-linux-gnu/ifort/


HDFINC = $(BASEDIR)/Linux/include/hdf/
HDFLIB = $(BASEDIR)/Linux/lib/

HDF5INC = $(BASEDIR)/Linux/include/hdf5/ 
HDF5LIB = $(BASEDIR)/Linux/lib/

HDF_INC = -I$(HDFINC) 
HDF_LIB = -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz -lm -L/usr/local/other/szip/2.1/noShared_intel-12.0.1.107/lib/ -lsz
HDF5_INC = -I$(HDF5INC)
HDF5_LIB = -L$(HDF5LIB) -lhdf5 -lhdf5_fortran -lhdf5


MYINC = $(HDF_INC)

PATH_OD_SHARED = 
PATH_CT_SHARED = 

CC = gcc 
CFLAGS = -O3 -funsigned-char -DMACINTOSH -DUSE_GDAS -c  

F77 = ifort
F77FLAGS = -O3  -c  
F90 = ifort 
F90FLAGS = -O3  
ADD_F90FLAGS = -c

OBJECT_LIST = modis_numerical_module.o \
			  general_array_io.o \
			  geos5_io_module.o \
			  global_model_grids.o \
			  modis_albedo.o \
			  emissivity.o \
			  read_albedo.o


all : read_albedo.exe

read_albedo.exe : $(OBJECT_LIST) 
	$(F90) $(F90FLAGS)  -o read_albedo.exe $(OBJECT_LIST) $(HDF_LIB) $(HDF5_LIB)
	cp read_albedo.exe ../MPI/

global_model_grids.o: $(PATH_OD_SHARED)global_model_grids.f90
	$(F90) $(F90FLAGS) $(ADD_F90FLAGS)  $(PATH_OD_SHARED)global_model_grids.f90
	
general_array_io.o : general_array_io.f90
	$(F90)  $(F90FLAGS)  $(HDF_INC) general_array_io.f90 -c

geos5_io_module.o : geos5_io_module.f90
	$(F90)  $(F90FLAGS)  $(HDF5_INC) geos5_io_module.f90 -c

modis_numerical_module.o : modis_numerical_module.f90
	$(F90) $(F90FLAGS) modis_numerical_module.f90 -c 
	
modis_albedo.o : modis_albedo.f90
	$(F90)  $(F90FLAGS) $(HDF_INC) modis_albedo.f90 -c

emissivity.o : emissivity.f90
	$(F90)  $(F90FLAGS)  $(HDF_INC) emissivity.f90 -c
	
read_albedo.o : read_albedo.f90
	$(F90)  $(F90FLAGS)  $(HDF_INC) $(HDF5_INC) read_albedo.f90 -c
	
