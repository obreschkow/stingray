# compiler
FC = gfortran

# compiler flags (try dfferent ones, if the default line does not work)
FCFLAGS = -g -O0 -fbounds-check -fwhole-file -ffpe-trap=invalid,zero,overflow -Wall -Wunused -Wuninitialized -Wsurprising -Wconversion
#FCFLAGS = -O3 -fopenmp

# List of executables to be built within the package
PROGRAMS = stingray

# "make" builds all
all: $(PROGRAMS)

stingray.o:    module_constants.o \
               module_types.o \
               module_system.o \
               module_parameters.o \
               module_cosmology.o \
               module_user.o \
               module_geometry.o \
               module_cone_intrinsic.o \
               module_cone_apparent.o \
               module_developer.o
stingray: 	   module_constants.o \
               module_types.o \
               module_system.o \
               module_parameters.o \
               module_cosmology.o \
               module_user.o \
               module_geometry.o \
               module_cone_intrinsic.o \
               module_cone_apparent.o \
               module_developer.o

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