include ../make.inc/Makefile.${PLAT}

DIRS =  $(PROJECT_BASE)/ocean/src\

EXE  = test
OBJS = test.o
OBJLIBS = lib$(LIB_OCEAN).a

# Header dependencies, especially important for templated stuff	
HEADERDEPS = \
	../newton/src/Newton.H \
	../newton/src/NewtonDecl.H \
	../thetastepper/src/ThetaStepper.H \
	../thetastepper/src/ThetaStepperDecl.H \
	../ocean/src/GlobalDefinitions.H \
	../ocean/src/Ocean.H \

INTLIBS = \
	-l$(LIB_OCEAN)\

all: $(EXE)

$(EXE):  $(OBJLIBS) $(OBJS)
	$(CXX) $(OBJS) -o $(EXE) $(LIBRARY_DIRS) $(INTLIBS) $(LIBS)

lib$(LIB_OCEAN).a: force_look
	cd $(PROJECT_BASE)/$(LIB_OCEAN)/src; $(MAKE) $(MFLAGS)

%.o: %.C $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

clean:
	-rm -f $(EXE) $(OBJS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look:	
	true
