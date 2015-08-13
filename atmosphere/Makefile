include ../make.inc/Makefile.${PLAT}

libbasename = $(LIB_ATMOS)
libname     = lib$(libbasename).a

obj= Atmosphere.o

HEADERDEPS = \
		../multiarray/MultiArray.H\
		../supervector/SuperVector.H\
		../ocean/GlobalDefinitions.H

all: lib

lib: $(libname)

$(libname): ${obj}
	ar vrcs $(libname) $(obj)

Atmosphere.o: Atmosphere.C Atmosphere.H
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

%.o: %.C %.H $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
%.o: %.C
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

clean:
	-rm -f *.o $(libname) *~ *#
