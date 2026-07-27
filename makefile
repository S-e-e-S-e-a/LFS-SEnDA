VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h  

#KHOME = /vol-th/software/lapack-3.3.1

#CF90 = /THL8/software/intel2018.4/impi/2018.4.274/bin64/mpiifort  #linrp 20200529
#CF90 = /public/software/mpi/intelmpi/2017.4.239/intel64/bin/mpiifort
CF90 = mpif90

CF77 = $(CF90)
LD = $(CF90)

#Paralellization opts
PARO =

#Size defaults
SIZEO = -r8 

#Arch opts
ARCHO= 

#Optimalization opts
OPTO= -O3 

#Inlining opts
INLO=

# Diverse options
DIVO= 

# Flags for debugging -- slow, gives location of SIGTRAP
DEBUG_FLAGS = 
DEBUG_FLAGS = -g

F77FLG = 
F90FLG = 

#NETCDFDIR  := /THL8/software/netcdf/4.4-icc16
#NETCDFDIR  := /public/software/mathlib/netcdf/4.4.1/intel_parallel
NETCDFDIR  := /public/software/mathlib/libs-intel/netcdf/4.4.1
INCLUDEDIR    := -I$(NETCDFDIR)/include
#FFLAGS    = -g -convert big_endian -assume byterecl -O3 -vec_report0  -shared-intel -mcmodel=large $(INCLUDEDIR)
FFLAGS    = -g -convert big_endian -assume byterecl -O3 -shared-intel -mcmodel=large $(INCLUDEDIR)
# linrp 20200529 remove vec_report0
# lyw remove vec_report0
FFLAGS_OE = $(FFLAGS)
LINKFLAGS =  $(FFLAGS)  

CPPARCH = 
CPPMODEL  =  -DPACIFIC 
CPPFLAGS = -P -traditional  $(CPPARCH) $(CPPMODEL) 
#LIBS =  -L$(NETCDFDIR)/lib/ -lnetcdf -L/THL8/software/lapack/3.8.0-icc16/lib64 -llapack -L/THL8/software/blas/3.7.0-icc16 -lblas  #linrp 20200529
#LIBS =  -L$(NETCDFDIR)/lib/ -lnetcdf -L/public/software/mathlib/lapack/3.8.0/intel/lib -llapack -L/public/software/compiler/intel/intel-compiler-2017.5.239/mkl/lib/intel64 -lblas
#LIBS =  -L$(NETCDFDIR)/lib/ -lnetcdf -L/public/software/mathlib/lapack/3.8.0/intel/lib -llapack -L/public/software/compiler/intel/intel-compiler-2017.5.239/mkl/lib/intel64 -lblas
#LIBS =  -L$(NETCDFDIR)/lib/ -lnetcdf -L /data/yuzp/LicomFS/assim/sla/Analysis -llapack -L /data/yuzp/LicomFS/assim/sla/Analysis -lblas
LIBS =  -L$(NETCDFDIR)/lib/ -lnetcdf -L ./ -llapack -L ./ -lblas
CPP = /usr/bin/cpp

# Rules for running cpp and updating files in TMP directory
.H.h:
	rm -f ./TMP/$*.h
	cat MODEL.CPP $*.H | $(CPP) $(CPPFLAGS) | ../sh/st.sh > ./TMP/$*.h


.F90.o:
	rm -f ./TMP/$*.f90
	cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS)  |../sh/st.sh > ./TMP/$*.f90
	cd ./TMP ; $(CF90) -c $(FFLAGS) $(F90FLG) -o $*.o $*.f90  

.F.o:
	rm -f ./TMP/$*.f
	cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) | ../sh/st.sh > ./TMP/$*.f
	cd ./TMP ; $(CF77) -c $(FFLAGS) $(F77FLG) -o $*.o $*.f  


TARGET =  enoi_zyz
TARGET1 = test_Lapack
TARGET2 = p_test
TARGET3 = prep_obs


include ../sh/source_zyz.files
#include omp_exception.files

INC2 =$(INC1:.H=.h)
FILES =$(F90FILES) $(F77FILES) $(MODULES)
FFILES =$(F90FILES:.F90=.f90) $(F77FILES:.F=.f) $(MODULES:.F90=.f90)
OBJECTS = $(F90FILES:.F90=.o) $(F77FILES:.F=.o) 
OMOD = $(MODULES:.F90=.o) $(MODULES77:.F=.o)

OMP_EXCEPTION_F77 =
OMP_EXCEPTION_F90 =
OMP_EXCEPTION_OBJ77=$(OMP_EXCEPTION_F77:.f=.o)
OMP_EXCEPTION_OBJ90=$(OMP_EXCEPTION_F90:.f90=.o)


all: $(TARGET) 


$(TARGET): $(INC2) $(OMOD) $(OBJECTS)  $(OMP_EXCEPTION_OBJ77) $(OMP_EXCEPTION_OBJ90)
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OMOD) $(OBJECTS) $(LIBS) 



# Rules for generating "omp exception" files
#$(OMP_EXCEPTION_F90) $(OMP_EXCEPTION_OBJ90):
#	rm -f ./TMP/$*.f90
#	cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS)  | ./st.sh > ./TMP/$*.f90 
#	cd ./TMP ; $(CF90) -c $(FFLAGS_OE) $(F90FLG) -o $*.o $*.f90

#$(OMP_EXCEPTION_F77) $(OMP_EXCEPTION_OBJ77):
#	rm -f ./TMP/$*.f
#	cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) | ./st.sh > ./TMP/$*.f
#	cd ./TMP ; $(CF77) -c $(FFLAGS_OE) $(F77FLG) -o $*.o $*.f

#################################################################################
OBJECTS1 =  p_test_Lapack.o
$(TARGET1): $(OBJECTS1) $(OMOD) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET1) $(OBJECTS1) $(LIBS) 

#################################################################################
OBJECTS2 =  p_test.o
$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS)


#################################################################################
OBJECTS3 =  mod_measurement.o mod_dimensions.o mod_angles.o mod_grid.o \
mod_modstate.o m_modstate_point.o m_modsubstate_point.o mod_meanssh.o \
.o m_read_mean_ssh.o p_prep_obs.o m_spherdist.o  m_prep_4_ENOI.o

$(TARGET3): $(OBJECTS3) $(OMOD)
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3) $(OBJECTS3) $(LIBS)

#################################################################################


clean:
	cd ./TMP ; rm *.f  *.o *.f90 *.h *.mod ../$(TARGET)


new: source depend

source:
	../sh/mksource.sh > ../sh/source_zyz.files

depend:
	../sh/mkdepend.pl | sort -u > ../sh/depends.file

omp_exceptions:
	./omp_except $(FFLAGS_OE) $(F90FLG)

include ../sh/depends.file

