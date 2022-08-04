# Makefile to build mpi class examples
# -1-- macros
#FC = $(F90)
FC = ifort -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

#--
# increasing the stacksize tio see if problem could be fixed
#--
MYFLAGS = -heap-arrays

#--
#trying to isolate the fault with a stack traceback
#---
OPT = -O3 -g -traceback -check all -fp-stack-check
#OPT = -O3 -g -traceback -fp-stack-check
#OPT = -O3 -fp-stack-check
OPT = -O3 -ipo -check bounds
#OPT = -O3 -ipo

#----
F90CFLAGS = $(MYFLAGS) $(OPT)
#F90CFLAGS = $(OPT)
F90LFLAGS = $(OPT)

EXE_SH = TiktakGlobalSearch

#----------
# -2- objects
#     List files from least USE to most USE
objects = nrtype.o stateControl.o genericParams.o utilities.o simplex.o objective.o minimize.o TiktakGlobalSearch.o

#--------------------------------------------
# -3- Define how to create the exe file
#--------------------------------------------
all: $(EXE_SH)

$(EXE_SH) : $(objects)
	@echo ""
	@echo "   LINKING to Create $@  ..."
	$(FC) $(F90LFLAGS) -o $@ $(objects)
	@echo ""
	@echo " -- Executable $@ created --"
	@echo ""

#--------------------------------------------
# -4- Define how create the .o files required
#     by the .exe file
#
#--------------------------------------------
nrtype.o: nrtype.f90
	$(FC) -c $(F90CFLAGS) $<

stateControl.o: stateControl.f90
	$(FC) -c $(F90CFLAGS) $<

genericParams.o: genericParams.f90
	$(FC) -c $(F90CFLAGS) $<

utilities.o: utilities.f90
	$(FC) -c $(F90CFLAGS) $<

simplex.o: simplex.f90
	$(FC) -c $(F90CFLAGS) $<

objective.o: objective.f90
	$(FC) -c $(F90CFLAGS) $<

minimize.o: minimize.f90
	$(FC) -c $(F90CFLAGS) $<

TiktakGlobalSearch.o:  TiktakGlobalSearch.f90
	$(FC) -c $(F90CFLAGS) $<

#---------------------------------------------
# -5- Define how to clean up the .o and .exe files
#     if you wish to build .exe from scratch
#---------------------------------------------
clean:
	 @echo "        Cleaning UP executable and object files ..."
	 \rm -f *.o *.mod $(objects) $(EXE_SH)
