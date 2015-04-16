include ../../make.inc/Makefile.${PLAT}

libbasename = $(LIB_OCEAN)
libname     = lib$(libbasename).a

obj= assemble.o atm.o monthly.o boundary.o forcing.o grid.o inout.o\
     itplbv.o levitus.o mat.o matetc.o lev.o mix.o res.o usr.o par.o\
     global.o thcm_utils.o scaling.o\
	 mix_imp.o mix_sup.o spf.o topo.o usrc.o\
	 Ocean.o\
     THCM.o OceanGrid.o\
	 TRIOS_Domain.o TRIOS_BlockPreconditioner.o TRIOS_Saddlepoint.o\
	 TRIOS_SolverFactory.o TRIOS_Static.o\
	 Utils.o

all: lib

test: $(libname) test.o 
	$(CXX) test.o  $(CXX_FLAGS) -o test $(LIBRARY_DIRS) -L. -l$(libbasename) $(LIBS)

test2: test2.o 
	$(CXX) test2.o  $(CXX_FLAGS) -o test2 $(LIBRARY_DIRS) $(LIBS)

test3: test3.o 
	$(CXX) test3.o  $(CXX_FLAGS) -o test3 $(LIBRARY_DIRS) $(LIBS)

test4: $(libname) test4.o 
	$(CXX) test4.o  $(CXX_FLAGS) -o test4 $(LIBRARY_DIRS) -L. -l$(libbasename) $(LIBS)

lib: $(libname)

$(libname): ${obj}
	ar vrcs $(libname) $(obj)

%.o: %.C %.H
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
%.o: %.C
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
%.o: %.F90 usr.o
	$(FC) -c $(FC_FLAGS) $<
%.o: %.f usr.o
	$(F77) -c $(F77_FLAGS) $<

# dependencies

par.o: par.F90
	$(FC) -c $(FC_FLAGS) $<

usr.o: usr.F90 par.o
	$(FC) -c $(FC_FLAGS) $<

global.o:      lev.o

assemble.o:    mat.o matetc.o atm.o

thcm_utils.o:  mix.o

boundary.o:    mat.o matetc.o atm.o

forcing.o:     inout.o itplbv.o matetc.o atm.o

inout.o:       global.o

levitus.o:     inout.o lev.o

scaling.o:     mat.o

usrc.o:        usrc.F90 mat.o atm.o res.o mix.o scaling.o monthly.o

matetc.o:      mat.o

mix_imp.o:     mat.o matetc.o mix_sup.o mix.o

mix_sup.o:     mix_sup.F90

spf.o:         atm.o

topo.o:        itplbv.o matetc.o global.o

clean:
	-rm -f *.o *.mod $(libname) *~ fort.* dump *.txt *#
