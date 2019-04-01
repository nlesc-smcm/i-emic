#                        Rules implementing Actions:
#                        ===========================
#
MRILU_ARCH=$(NUMBASE)/$(ARCH)
MRILU_BIN=$(MRILU_ARCH)/bin
MRILU_LIB=$(MRILU_ARCH)/lib
MRILU_MOD=$(MRILU_ARCH)/mod 
MRILU_OBJ=$(MRILU_ARCH)/obj
MRILU_DOC=$(NUMBASE)/doc
# Creating all directories

.PHONY: dirs
#
dirs: $(NUMBASE) $(MRILU_ARCH) $(MRILU_BIN) $(MRILU_LIB) $(MRILU_MOD) $(MRILU_OBJ) $(MRILU_DOC)
#
# Rules to create a specific directory
#
$(NUMBASE):
	mkdir $(NUMBASE)
	-chmod go+rx $(NUMBASE)
#
$(MRILU_ARCH):
	mkdir $(MRILU_ARCH)
	-chmod go+rx $(MRILU_ARCH)
#
$(MRILU_BIN):
	mkdir $(MRILU_BIN)
	-chmod go+rx $(MRILU_BIN)
#
$(MRILU_LIB):
	mkdir $(MRILU_LIB)
	-chmod go+rx $(MRILU_LIB)
#
$(MRILU_MOD):
	mkdir $(MRILU_MOD)
	-chmod go+rx $(MRILU_MOD)
#
$(MRILU_OBJ):
	mkdir $(MRILU_OBJ)
	-chmod go+rx $(MRILU_OBJ)
#
$(MRILU_DOC):
	mkdir $(MRILU_DOC)
	-chmod go+rx $(MRILU_DOC)

# Disable parallel make
.NOTPARALLEL:

#
# Default goal: compile all modules
#
.PHONY: default_all
#
default_all:\
            $(addsuffix .mod, $(addprefix m_, $(modules) ) )\
            $(addsuffix .o, $(progs) )\
            $(addsuffix .o, $(objects) )
#
# Installation
#
.PHONY: default_install
#
default_install: dirs install_others install_doc all install_obj install_prg
#
# Making the documentation
#
.PHONY: default_install_doc
#
default_install_doc:\
             $(addprefix $(NUMBASE)/doc/, $(addsuffix .txt, $(modules) ) ) \
             $(addprefix $(NUMBASE)/doc/, $(addsuffix .txt, $(progs) ) ) \
             $(addprefix $(NUMBASE)/doc/, $(addsuffix .txt, $(objects) ) )
#
#
# Adding the object files to the library
#
.PHONY: default_install_obj
#
default_install_obj:\
	$(MRILU_LIB)/$(thislib).a($(addsuffix .o, $(modules)) $(addsuffix .o, $(objects)) )
#
.PHONY: default_install_prg
#
install_prg:\
             $(addprefix $(MRILU_BIN)/, $(progs) )\
#


.PHONY: default_uninstall
#
default_uninstall: default_clean default_clean_doc default_clean_prog

.PHONY: default_clean
#
default_clean: default_clean_bak default_clean_obj
#
.PHONY: default_clean_bak
#
default_clean_bak: 
	$(RM) -r *~
#	$(RM) -r .nfs*
	$(RM) -r \#*\#
#
.PHONY: default_clean_obj
#
default_clean_obj: 
	$(RM) $(MRILU_OBJ)/*.o

.PHONY: default_clean_prog
#
default_clean_prog: 
	$(RM) $(addprefix $(MRILU_BIN)/, $(progs))
	$(RM) $(addprefix $(MRILU_MOD)/m_, $(addsuffix .mod , $(modules)))
	$(AR) -d $(MRILU_LIB)/$(thislib).a $(addsuffix .o, $(modules))
#
.PHONY: default_clean_doc
#
default_clean_doc:
	$(RM) $(addsuffix .txt, $(progs) )
	$(RM) $(addsuffix .txt, $(modules) )

#
#                        Pattern rules:
#                        ==============
#
# Rules to indicate the extra dependencies of the executables on the
# library routines:
#
#
####LDLIBES+=mkl_blas95 
###EXTRA= -L$(MKLPATH) -I$(MKLINCLUDE) -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -lpthread
$(MRILU_BIN)/%: %.o
	echo MKL=$(MKL)
	echo $(LDLIBDIR)
	echo $(LDLIBES)
	$(LD) $(LDFLAGS) $(MRILU_OBJ)/$*.o $(addprefix -L,$(LDLIBDIR)) $(addprefix -l,$(LDLIBES)) -o $@
#
#
#
#
#
# Pattern rule to update object file in library from the '.o' file:
#
(%.o): %.o
	$(AR) $(ARFLAGS) $@ $(MRILU_OBJ)/$*.o
#	echo archive $@
	-chmod u=rw,g=rw,o=r $@
#
#
# Pattern rule to produce a '.o' file from a '.c' file.
#
%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $*.c -o $(MRILU_OBJ)/$*.o
#
#%.o: m_%.mod
#
# Pattern rule to produce a '.o' file from a '.f90' file.
#
%.o: %.f90
	$(F90) $(FFLAGS) -c -I$(MRILU_MOD) $(MODFLAG)$(MRILU_MOD) $*.f90 -o $(MRILU_OBJ)/$*.o
#
# Pattern rule to produce a '.o' file from a '.f90.fpp' file.
#
m_%.mod: %.F90
	$(F90) $(FPPFLAGS) $(FFLAGS) -c -I$(MRILU_MOD)  $(MODFLAG)$(MRILU_MOD) $*.F90 -o $(MRILU_OBJ)/$*.o
#
# Pattern rule to produce a '.o' file from a '.f90.fpp' file.
#
%.o: %.F90
	$(F90) $(FPPFLAGS) $(FFLAGS) -c -I$(MRILU_MOD)  $(MODFLAG)$(MRILU_MOD) $*.F90 -o $(MRILU_OBJ)/$*.o

# Pattern rule to produce a '.o' file from a '.f90.fpp' file.
#
%.o: %.f90
	$(F90) $(FFLAGS) -c -I$(MRILU_MOD)  $(MODFLAG)$(MRILU_MOD) $*.F90 -o $(MRILU_OBJ)/$*.o
#
#
# Pattern rule to produce a '.o' file from a '.F' file.
#
%.o: %.F
	$(F77) $(FPPFLAGS) $(FFLAGS) -c -I$(MRILU_MOD)  $(MODFLAG)$(MRILU_MOD) $*.F -o $(MRILU_OBJ)/$*.o
#
#
# Pattern rule to produce a '.o' file from a '.f' file.
#
%.o: %.f
	$(F77) $(FFLAGS) -c -I$(MRILU_MOD)  $(MODFLAG)$(MRILU_MOD) $*.f -o $(MRILU_OBJ)/$*.o
#
#
# Pattern rule to extract the documentation file and place it in the
# current directory:
#
%.txt: %.c
	@echo "NO documentation for .c files"
#
%.txt: %.F
	$(MRILU_BIN)/getdoc $< $@
#
%.txt: %.f
	$(MRILU_BIN)/getdoc $< $@
#
#
# Pattern rule to extract the documentation file and place it in the
# directory  $(NUMBASE)/doc/:
#
$(NUMBASE)/doc/%.txt: %.c
	@echo "NO documentation for .c files"
#
$(NUMBASE)/doc/%.txt: %.F90
	$(RM) $@
	$(MRILU_BIN)/getdoc $< $@
	-chmod u=rw,g=rw,o=r $@
#
$(NUMBASE)/doc/%.txt: %.f90
	$(RM) $@
	$(MRILU_BIN)/getdoc $< $@
	-chmod u=rw,g=rw,o=r $@
$(NUMBASE)/doc/%.txt: %.F
	$(RM) $@
	$(MRILU_BIN)/getdoc $< $@
	-chmod u=rw,g=rw,o=r $@
#
$(NUMBASE)/doc/%.txt: %.f
	$(RM) $@
	$(MRILU_BIN)/getdoc $< $@
	-chmod u=rw,g=rw,o=r $@
#
# installing 
#
#$(MRILU_LIB)/libmkl_blas95.a:
#	cd $(MKL_PATH)/interfaces/blas95; $(MAKE) lib PLAT=lnx32
#	mv -f $(MKL_PATH)/interfaces/blas95/libmkl_blas95.a $@
#	-chmod u=rw,g=rw,o=r $@
#
#$(MRILU_MOD)/mkl95_blas.mod:
#	cd $(MKL_PATH)/interfaces/blas95; $(MAKE) lib PLAT=lnx32
#	mv -f $(MKL_PATH)/interfaces/blas95/mkl95_blas.mod $@
#	-chmod u=rw,g=rw,o=r $@
#
#$(MRILU_MOD)/mkl95_precision.mod:
#	cd $(MKL_PATH)/interfaces/blas95; $(MAKE) lib PLAT=lnx32
#	mv -f $(MKL_PATH)/interfaces/blas95/mkl95_precision.mod $@
#	-chmod u=rw,g=rw,o=r $@
#












