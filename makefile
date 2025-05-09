# Call as make [sam=shark,...] [mode=standard,dev]

# optional arguments:
# sam = galaxy formation model; each model requires custom modules "module_user_routines_[sam].f90" and "module_user_selection_[sam].f90"
# mode = compilation mode; allowed modes are 'default' and 'dev'

# NB: If the HDF5 library cannot be found by the compiler, try to set the environment variables
#     HDF5_INC and HDF5_LIB to the paths needed for the compiler flags -I and -L

ifndef sam
   sam = shark
endif

ifndef mode
   mode = default
endif

# HDF5 library flags
LFLAGS = -lhdf5_fortran -lhdf5
ifneq ($(HDF5_LIB),)
    LFLAGS += -L$(HDF5_LIB)
endif
ifneq ($(HDF5_INC),)
    LFLAGS += -I$(HDF5_INC)
endif

# Other compiler flags
ifeq ($(mode),default)
   CFLAGS = -O3 -fopenmp -ffree-line-length-0
else ifeq ($(mode),dev)
   CFLAGS = -O0 -g -fbounds-check -fwhole-file -ffpe-trap=invalid,zero,overflow -Wall -Wunused -Wuninitialized -Wsurprising -Wconversion
else
   $(error ERROR: unknown mode: '${mode}')
endif

# concatenate flags
FCFLAGS =  $(CFLAGS) $(LFLAGS)

# Compiler
FC = gfortran

# user info
$(info Compilation options:)
$(info + Computing system = '${system}'.)
$(info + Compiling mode = '${mode}'.)
$(info + Galaxy formation model = '${sam}'.)

# List of executables to be built within the package
PROGRAMS = stingray

# "make" builds all
all: $(PROGRAMS)

stingray.o:    shared/shared_module_core.o \
               shared/shared_module_arguments.o \
               shared/shared_module_parameters.o \
               shared/shared_module_hdf5.o \
               shared/shared_module_cosmology.o \
               shared/shared_module_maths.o \
               shared/shared_module_vectors.o \
               shared/shared_module_constants.o \
               shared/shared_module_sort.o \
               module_global.o \
               module_parameters.o \
               module_conversion.o \
               module_emission_lines.o \
               module_user_routines_$(sam).o \
               module_selection_tools.o \
               module_user_selection_$(sam).o \
               module_tiling.o \
               module_sky.o
               
stingray: 	   shared_module_core.o \
               shared_module_arguments.o \
               shared_module_parameters.o \
               shared_module_hdf5.o \
               shared_module_cosmology.o \
               shared_module_maths.o \
               shared_module_vectors.o \
               shared_module_constants.o \
               shared_module_sort.o \
               module_global.o \
               module_parameters.o \
               module_conversion.o \
               module_emission_lines.o \
               module_user_routines_$(sam).o \
               module_selection_tools.o \
               module_user_selection_$(sam).o \
               module_tiling.o \
               module_sky.o

# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) -o $@ $^ $(FCFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
	
%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

%.o: %.f03
	$(FC) $(FCFLAGS) -c $<
	
%.o: %.f08
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD
	rm -f *~ $(PROGRAMS)
	rm -f fort.*
	rm -rf *.dSYM
	rm -rf example
