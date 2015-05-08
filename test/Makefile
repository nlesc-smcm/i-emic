include ../make.inc/Makefile.${PLAT}

DIRS =  $(PROJECT_BASE)/$(LIB_OCEAN)/

EXE  = test timestepping
OBJS = test.o timestepping.o
OBJLIBS = lib$(LIB_OCEAN).a

# Header dependencies, especially important for templated stuff	
HEADERDEPS = \
	../newton/Newton.H \
	../newton/NewtonDecl.H \
	../vector/Vector.H\
	../thetastepper/ThetaStepper.H \
	../thetastepper/ThetaStepperDecl.H \
	../continuation/Continuation.H \
	../continuation/ContinuationDecl.H \
	../ocean/GlobalDefinitions.H \
	../ocean/Ocean.H \

INTLIBS = \
	-l$(LIB_OCEAN)\

CXX_DB = -DDEBUGGING

all: test

test: $(OBJLIBS) test.o
	$(CXX) test.o -o test $(LIBRARY_DIRS) $(INTLIBS) $(LIBS)

timestepping: $(OBJLIBS) timestepping.o
	$(CXX) timestepping.o -o timestepping $(LIBRARY_DIRS) $(INTLIBS) $(LIBS)

lib$(LIB_OCEAN).a: force_look
	cd $(PROJECT_BASE)/$(LIB_OCEAN); $(MAKE) $(MFLAGS)

%.o: %.C $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $(CXX_DB) $<

clean:
	-rm -f $(EXE) $(OBJS) $(OBJLIBS) fort.* *.txt *.eps
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look:	
	true
