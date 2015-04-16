include ../make.inc/Makefile.${PLAT}

DIRS =  $(PROJECT_BASE)/$(LIB_OCEAN)/src\
        $(PROJECT_BASE)/$(LIB_NEWTON)/src	

EXE  = test
OBJS = test.o
OBJLIBS = \
	lib$(LIB_OCEAN).a	\
    lib$(LIB_NEWTON).a	\

# Header dependencies, especially important for templated stuff	
HEADERDEPS = \
	../newton/src/Newton.H\
	../ocean/src/GlobalDefinitions.H\
	../ocean/src/Ocean.H\

INTLIBS = \
	-l$(LIB_OCEAN)	\
	-l$(LIB_NEWTON)	\

all: $(EXE)

$(EXE): $(OBJS) $(OBJLIBS) 
	$(CXX) $(OBJS) -o $(EXE) $(LIBRARY_DIRS) $(INTLIBS) $(LIBS)

lib$(LIB_OCEAN).a: force_look
	cd $(PROJECT_BASE)/$(LIB_OCEAN)/src; $(MAKE) $(MFLAGS)

lib$(LIB_NEWTON).a: force_look
	cd $(PROJECT_BASE)/$(LIB_NEWTON)/src; $(MAKE) $(MFLAGS)

$(HEADER):

%.o: %.C $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

clean:
	-rm -f $(EXE) $(OBJS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look:	
	true
