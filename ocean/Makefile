include ../make.inc/Makefile.${PLAT}

libbasename = $(LIB_OCEAN)
libname     = lib$(libbasename).a

obj= assemble.o atm.o monthly.o boundary.o forcing.o grid.o inout.o\
     itplbv.o levitus.o mat.o matetc.o lev.o mix.o res.o usr.o par.o\
     global.o thcm_utils.o scaling.o\
	 mix_imp.o mix_sup.o spf.o topo.o usrc.o\
	 inserts.o\
	 Ocean.o OceanTheta.o OceanCont.o\
     THCM.o OceanGrid.o\
	 TRIOS_Domain.o TRIOS_BlockPreconditioner.o TRIOS_Saddlepoint.o\
	 TRIOS_SolverFactory.o TRIOS_Static.o\
	 Utils.o

HEADERDEPS = \
	$(PROJECT_BASE)/vector/Vector.H\
	$(PROJECT_BASE)/ocean/GlobalDefinitions.H

all: lib

lib: $(libname)

$(libname): ${obj}
	ar vrcs $(libname) $(obj)

Ocean.o: Ocean.C Ocean.H $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
OceanTheta.o: OceanTheta.C OceanTheta.H $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
OceanCont.o: OceanCont.C OceanCont.H $(HEADERDEPS)
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<

%.o: %.C %.H
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
%.o: %.C
	$(CXX) -c $(INCLUDE_DIRS) $(CXX_FLAGS) $<
%.o: %.F90 usr.o
	$(FC) -c $(FC_FLAGS) $<
%.o: %.f usr.o
	$(F77) -c $(F77_FLAGS) $<

# specific dependencies

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
