include ./make.inc/Makefile.${PLAT}

all:
	cd utils/src;   make -j $(JOBS)
	cd ocean/src;   make lib -j $(JOBS)
	cd newton/src;  make -j $(JOBS)
