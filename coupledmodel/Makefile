include ../make.inc/Makefile.${PLAT}

libbasename = $(LIB_COUPLED)
libname     = lib$(libbasename).a

obj= CoupledModel.o

HEADERDEPS = \
	../supervector/SuperVector.H\
    ../atmosphere/Atmosphere.H\
    ../idrsolver/IDRSolver.H\
    ../gmressolver/GMRESSolver.H\
    ../ocean/Ocean.H\

all: lib

lib: $(libname) 

$(libname): ${obj} 
	ar vrcs $(libname) $(obj)

CoupledModel.o: CoupledModel.C CoupledModel.H $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

%.o: %.C %.H 
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
%.o: %.C
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

clean:
	-rm -f *.o $(libname) *~ *#
