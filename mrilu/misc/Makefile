#                                                            2002-10-07
#
#                        Description:
#                        ============
#
# Makefile for the subroutines in  misc/.
# Make can be used to:
# . produce or update the object files in $(MRILU_LIB)/libgen.a
# . produce or update the documentation files in  $(NUMBASE)/doc/
#
#
# Usage:
#
#   make [ OPT = [ debug | optim | prof ] ]
#        [ CC = <name C-Compiler> ] [ F77 = <name Fortran 77-Compiler ]
#        <action>
#
#   The contents of the C and Fortran 77 source files depend on the value of
#   the OPTion  variable definition, OPT.
#   The preprocessor(s) and compiler(s) depend on the value of the
#   the options CC or F77
#   The compilation options depend on the value of the  OPTion  variable
#   definition, OPT, and the compilers to be used.
#
#   For the different definitions we have:
#   OPT=  or omitted
#     The compiler optimisation option '-O' is used to produce the object
#     files.  The range checking and profiling are disabled.
#   OPT=debug
#     The preprocessor symbol DEBUG is defined.
#     The compiler optimisation option is disabled and the options for
#     range checking and profiling are enabled.  Extra trace code is
#     inserted.
#   OPT=optim
#     The preprocessor symbol OPTIM is defined.
#     The compiler optimisation options are enabled to produce the object
#     files.  Range checking and profiling are disabled.
#   OPT=prof
#     The profiling is enabled.  The compiler optimisation option is enabled,
#     if possible, and the range checking is disabled.
#   CC=<name C-compiler>
#     The name of C compiler is <name C-compiler>.
#   F77=<name Fortran 77-compiler>
#     The name of Fortran 77 compiler is <name Fortran 77-compiler>.
#
#   The following actions are defined:
#   all  (the default action)
#     produce the object files and the executable test programs.
#   clean
#     remove all redundant files.
#   clean_doc
#     remove all the document files (extension '.txt') from the current
#     directory.
#   clean_ftn
#     remove each '.f', with an existing pre-Fortran '.F' file, from the
#     current directory.
#   clean_obj
#     remove all the object files from the current directory.
#   clean_test
#     remove all the executable test programs with the corresponding output
#     files.
#   doc
#     extract the documentation from the source files.
#   fortran_src
#     produce Fortran source files, with extension '.f', from the pre-Fortran
#     source files with extension '.F'.
#   install
#     update the document files in    $(NUMBASE)/doc/,
#            the object files in      $(MRILU_LIB)/libgen.a
#   install_doc
#     place the document files of the changed source files in the directory
#     $(NUMBASE)/doc/
#   install_obj
#     place the object files of the changed source files in the library
#     $(MRILU_LIB)/libgen.a
#   test
#     generate the testprograms and produce the test output files.
#
#
# N.B.
# The make file uses the following environment variables:
#   NUMBASE  The base of the subdirectories in which binaries, libraries
#            etc. reside.

#
#                        Definitions:
#                        ============
#
# include userdefined constants
include ../makefile.inc
#
#
# Define the base names of the source files of the subroutines involved:
#
modules =  dgeco   dgedi   dgefa   dgesl  dperv iperv csr2cscvv

thislib = libmisc
#
#
ifdef LIBATLAS
FFLAGS  += -DWITH_ATLAS
else
modules += dgetrf dgetrs dgetf2 dlaswp dtrsm dswap
dswap   =  m_dswap.mod
endif
#
 Define the base names of the source files of the test programs involved:
#
#
# Definition of search path variable:
#
VPATH = $(MRILU_LIB):$(MRILU_MOD):$(MRILU_OBJ)
#
#
#
#
# Delete target in case of errors during construction of target
#
.DELETE_ON_ERROR:
#
#
.SUFFIXES:                          # Delete the default suffixes
.SUFFIXES: .o .c .f .F .txt         # Define the used suffix list
#
#
#
#                    Pattern rules:
#                    ==============
#
#
include $(SOURCEDIR)/mrilu/scripts/rules.mk
#
# Pattern rules to extract the documentation files and place them in the
# current directory:
#
%.txt: %.c
	@echo "NO documentation for .c files"
%.txt: %.F
	$(NUMBASE)/bin/getdoc $< $@
%.txt: %.f
	$(NUMBASE)/bin/getdoc $< $@
#
#
# Pattern rules to extract the documentation files and place them in
# $(NUMBASE)/doc/:
#
$(NUMBASE)/doc/%.txt: %.c
	@echo "NO documentation for .c files"
$(NUMBASE)/doc/%.txt: %.F
	$(RM) $@
	$(NUMBASE)/bin/getdoc $< $@
	-chmod u=rw,g=rw,o=r $@
$(NUMBASE)/doc/%.txt: %.f
	$(RM) $@
	$(NUMBASE)/bin/getdoc $< $@
	-chmod u=rw,g=rw,o=r $@
#
#
# Pattern rules to run a test program and generate an output file.
#
%.out: %.F
	$* > $@
%.out: %.f
	$* > $@
#
#
#                        Rules implementing Actions:
#                        ===========================
#
#
# Default goal: compile all modules
#
.PHONY: all
#
all: default_all
#
# Installation
#
.PHONY: install
#
install: default_install
#
# Making the documentation
#
.PHONY: install_doc
#
install_doc: default_install_doc
#
# Adding the object files to the library
#
.PHONY: install_obj
#
install_obj: default_install_obj
#
#


.PHONY: install_others
install_others:
#	cd ../getdoc;     $(MAKE) install
#
#
# Program  tstcopy
#
objtstcopy = dsprav.o icopy.o isprav.o tstcopy.o
tstcopy:  $(objtstcopy)
	$(FC) -o tstcopy $(LDFLAGS) $(objtstcopy)
#
#
#
.PHONY: uninstall
#
uninstall: default_uninstall

.PHONY: clean
#
clean: default_clean
#
.PHONY: clean_bak
#
clean_bak: default_clean_bak
#
.PHONY: clean_obj
#
default_clean_obj: clean_obj 

.PHONY: clean_prog
#
clean_prog: default_clean_prog

.PHONY: default_clean_doc
#
clean_doc: default_clean_doc

#
#
#
.PHONY: doc
#
doc: $(addsuffix .txt, $(modules) )
#
#
#                        Other Rules:
#                        ============
#
m_dgeco.mod:	m_dump.mod	m_dgefa.mod
m_dgedi.mod:	m_dump.mod	$(dswap)
m_dgetf2.mod:	$(dswap)
m_dgetrf.mod:	m_dump.mod	m_dtrsm.mod	m_dlaswp.mod	m_dgetf2.mod
m_dgetrs.mod:	m_dump.mod	m_dtrsm.mod	m_dlaswp.mod
m_dlaswap.mod:	$(dswap)
m_dperv.mod:	m_dump.mod
m_dtrsm.mod:	m_dump.mod
m_iperv.mod:	m_dump.mod
m_csr2cscvv.mod: m_wacsr.mod
#
#
# End of Makefile
