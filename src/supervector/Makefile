include ../make.inc/Makefile.${PLAT}

libbasename = $(LIB_SUPERVEC)
libname     = lib$(libbasename).a

obj= SuperVector.o

HEADERDEPS = \
		../ocean/GlobalDefinitions.H\

all: lib

lib: $(libname)

$(libname): ${obj}
	ar vrcs $(libname) $(obj)

SuperVector.o: SuperVector.C SuperVector.H
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

%.o: %.C %.H $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
%.o: %.C
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

clean:
	-rm -f *.o $(libname) *~ *#
