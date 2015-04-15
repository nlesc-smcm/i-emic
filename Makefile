
numprocs = 4

all:
	cd utils/src;   make -j $(numprocs)
	cd ocean/src;   make lib -j $(numprocs)
	cd newton/src;  make -j $(numprocs)
