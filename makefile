# Input model; each model requires a custom module "module_user_[SAM].f95"
Model = shark

# Compiler
FC = gfortran

# Compiler flags for testing and debugging
FCFLAGS = -g -O0 -I/usr/local/lib/hdf5/include -L/usr/local/lib/hdf5/lib -lhdf5_fortran -lhdf5 -fbounds-check -fwhole-file -ffpe-trap=invalid,zero,overflow -Wall -Wunused -Wuninitialized -Wsurprising -Wconversion

# Compiler flags for optimize exectusion
#FCFLAGS = -O3 -fopenmp -I/usr/local/lib/hdf5/include -L/usr/local/lib/hdf5/lib -lhdf5_fortran -lhdf5

# List of executables to be built within the package
PROGRAMS = stingray

# "make" builds all
all: $(PROGRAMS)

stingray.o:    module_types.o \
               module_constants.o \
               module_system.o \
               module_cosmology.o \
               module_conversion.o \
               module_hdf5_utilities.o \
               module_user_$(Model).o \
               module_parameters.o \
               module_geometry.o \
               module_cone_intrinsic.o \
               module_cone_apparent.o
               
stingray: 	   module_types.o \
               module_constants.o \
               module_system.o \
               module_cosmology.o \
               module_conversion.o \
               module_hdf5_utilities.o \
               module_user_$(Model).o \
               module_parameters.o \
               module_geometry.o \
               module_cone_intrinsic.o \
               module_cone_apparent.o

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