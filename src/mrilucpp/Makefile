include ../make.inc/Makefile.${PLAT}

MRILU_LIBDIR=-L${MRILUDIR}/lib
MRILU_MODDIR=-I${MRILUDIR}/mod

CXX_OBJ=Ifpack_MRILU.o
F_OBJ=mrilucpp.o

INCLUDE_DIRS+=${MRILU_MODDIR}

CXX_FLAGS+=-DHAVE_IFPACK_MRILU
FF_FLAGS+=-DHAVE_IFPACK_MRILU

all: lib

lib: libmrilucpp.a libifpack_mrilu.a

mrilucpp.o: mrilucpp.F90
	${FC} ${MRILU_MODDIR} ${FF_FLAGS} ${FC_FLAGS} -c $<

Ifpack_MRILU.o: Ifpack_MRILU.cpp Ifpack_MRILU.h
	${CXX} ${INCLUDE_DIRS} ${CXX_FLAGS} ${HEADERDEPS} -c $<

libmrilucpp.a: mrilucpp.o 
	ar vr libmrilucpp.a mrilucpp.o 

libifpack_mrilu.a: ${CXX_OBJ}
	ar vr libifpack_mrilu.a ${CXX_OBJ} 

clean:
	-rm -rf *.o *.mod *.a 
