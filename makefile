# Input model; each model requires a custom module "module_user_[SAM].f03"
Model = shark

# Compiler
FC = gfortran

# Compiler flags for testing and debugging
#FCFLAGS = -g -O0 -I/usr/local/lib/hdf5/include -L/usr/local/lib/hdf5/lib -lhdf5_fortran -lhdf5 -fbounds-check -fwhole-file -ffpe-trap=invalid,zero,overflow -Wall -Wunused -Wuninitialized -Wsurprising -Wconversion
#FCFLAGS = -g -O0 -I/usr/local/include -L/usr/local/lib -lhdf5_fortran -lhdf5 -fbounds-check -fwhole-file -ffpe-trap=invalid,zero,overflow -Wall -Wunused -Wuninitialized -Wsurprising -Wconversion

# Compiler flags for optimized execution
# FCFLAGS = -O3 -fopenmp -I/usr/local/lib/hdf5/include -L/usr/local/lib/hdf5/lib -lhdf5_fortran -lhdf5
FCFLAGS = -O3 -fopenmp -I/usr/local/include -L/usr/local/lib -lhdf5_fortran -lhdf5

# Compiler flags on hyades
FCFLAGS = -g -O3 -fopenmp -I/opt/bldr/local/storage/hdf5/1.10.2/include -L/opt/bldr/local/storage/hdf5/1.10.2/lib -lhdf5_fortran -lhdf5 -ffree-line-length-none

# List of executables to be built within the package
PROGRAMS = stingray

# "make" builds all
all: $(PROGRAMS)

stingray.o:    module_constants.o \
               module_system.o \
               module_types.o \
               module_io.o \
               module_linalg.o \
               module_cosmology.o \
               module_sort.o \
               module_conversion.o \
               module_hdf5.o \
               module_user_routines_$(Model).o \
               module_user_selections_$(Model).o \
               module_parameters.o \
               module_tiling.o \
               module_sky.o
               
stingray: 	   module_constants.o \
               module_system.o \
               module_types.o \
               module_io.o \
               module_linalg.o \
               module_sort.o \
               module_cosmology.o \
               module_conversion.o \
               module_hdf5.o \
               module_user_routines_$(Model).o \
               module_user_selections_$(Model).o \
               module_parameters.o \
               module_tiling.o \
               module_sky.o

# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
	
%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

%.o: %.f03
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD
	rm -f *~ $(PROGRAMS)
	rm -f fort.*
