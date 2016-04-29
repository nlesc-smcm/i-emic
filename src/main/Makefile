include ../make.inc/Makefile.${PLAT}

DIRS =  \
	$(PROJECT_BASE)/ocean/\
	$(PROJECT_BASE)/atmosphere/\
	$(PROJECT_BASE)/supervector/\
	$(PROJECT_BASE)/mrilucpp/\
	$(PROJECT_BASE)/coupledmodel/\

EXE     = run_coupled run_ocean test test_numjac
OBJS    = run_coupled.o run_ocean.o test.o test_numjac
OBJLIBS = lib$(LIB_OCEAN).a\
	  lib$(LIB_ATMOS).a\
	  lib$(LIB_COUPLED).a\
	  lib$(LIB_SUPERVEC).a\
      lib$(LIB_MRILUCPP).a\

# Header dependencies, especially important for templated stuff
HEADERDEPS = \
	../newton/Newton.H \
	../newton/NewtonDecl.H \
	../supervector/SuperVector.H\
	../thetastepper/ThetaStepper.H \
	../thetastepper/ThetaStepperDecl.H \
	../continuation/Continuation.H \
	../continuation/ContinuationDecl.H \
	../idrsolver/IDRSolver.H \
	../gmressolver/GMRESSolver.H \
	../ocean/GlobalDefinitions.H \
	../ocean/Ocean.H \
	../atmosphere/Atmosphere.H\
	../coupledmodel/CoupledModel.H\
	../main/NumericalJacobian.H\

# Coupled depends on ocean and atmos
INTLIBS = \
	-l$(LIB_COUPLED)\
	-l$(LIB_SUPERVEC)\
	-l$(LIB_OCEAN)\
	-l$(LIB_ATMOS)\

all: run_coupled

run_coupled: $(OBJLIBS)  run_coupled.o
	$(CXX) run_coupled.o -o run_coupled $(LIBRARY_DIRS) $(INTLIBS) $(LIBS)

run_ocean: $(OBJLIBS) run_ocean.o 
	$(CXX) run_ocean.o -o run_ocean $(LIBRARY_DIRS) $(INTLIBS) $(LIBS)

test: $(OBJLIBS) test.o
	$(CXX) test.o -o test $(LIBRARY_DIRS) $(INTLIBS) $(LIBS) -lgtest -lpthread

test_numjac: $(OBJLIBS) test_numjac.o
	$(CXX) test_numjac.o -o test_numjac $(LIBRARY_DIRS) $(INTLIBS) $(LIBS)

lib$(LIB_OCEAN).a: force_look
	cd $(PROJECT_BASE)/$(LIB_OCEAN); $(MAKE) $(MFLAGS)

lib$(LIB_ATMOS).a: force_look
	cd $(PROJECT_BASE)/$(LIB_ATMOS); $(MAKE) $(MFLAGS)

lib$(LIB_COUPLED).a: force_look
	cd $(PROJECT_BASE)/$(LIB_COUPLED); $(MAKE) $(MFLAGS)

lib$(LIB_SUPERVEC).a: force_look
	cd $(PROJECT_BASE)/$(LIB_SUPERVEC); $(MAKE) $(MFLAGS)

lib$(LIB_MRILUCPP).a: force_look
	cd $(PROJECT_BASE)/$(LIB_MRILUCPP); $(MAKE) $(MFLAGS)

%.o: %.C $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $(CXX_DB) $<

clean:
	-rm -f $(EXE) $(OBJS) $(OBJLIBS) fort.*  \
	*.eps *.co *.jco *.beg *.rl *.info *~ dump *#
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

include:
	echo $(INCLUDE_DIRS) 

force_look:
	true
