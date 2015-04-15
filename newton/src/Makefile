include ../../make.inc/Makefile.${PLAT}	

libbasename = $(LIB_NEWTON)
libname     = lib$(libbasename).a

obj = Newton.o

default: lib

lib: $(libname)

$(libname): ${obj}
	ar vrcs $(libname) $(obj)

%.o: %.C %.H
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

clean:
	rm *.o $(libname) *~
