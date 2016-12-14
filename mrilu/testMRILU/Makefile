#                                                            2003-02-25
#
#                        Description:
#                        ============
#
# Makefile for a main program, which uses the preconditioner software.
# 
# Make can be used to:
# . generate the program to test the subroutines in the libraries
#   $(MRILU_LIB)/libgen.a  and   $(MRILU_LIB)/libprecon.a
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
#   options CC or F77
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
#     range checking and profiling are enabled to produce the object
#     files.
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
#     generate the testprogram(s).
#   clean
#     remove all redundant files.
#   clean_obj
#     remove all the object files from the current directory.
#
#
# N.B.
# The make file uses the following environment variables:
#   NUMBASE  The base of the subdirectories in which binaries, libraries
#            etc. reside.

#
#
#                        Definitions:
#                        ============
# include userdefined constants
include ../makefile.inc
#
# Define the base names of the source files of the subroutines involved:
#
modules = arthura cmpsolprc eugen2d readsherman cmpsolprc2
#
#
# Define the base names of the source files of the programs involved:
#
progs = nonsym smallnonsym nschmat
thislib = libtestmrilu
#
#
#
# Load from libraries:
#
LD       = $(F90)
LDLIBDIR = $(MRILU_LIB)
LDLIBES  += testmrilu precon mtstor misc iosrc
#
ifeq "$(MKL)" "1"
LDLIBDIR += $(MKL_LIB)
LDLIBES  += mkl_blas95 mkl
mkl95_blas = libmkl_blas95.a
endif
#
ifdef $(LIBATLAS)
FFLAGS   += -DWITH_ATLAS
LDLIBDIR += $(LIBATLAS)
LDLIBES  += lapack cblas f77blas atlas g2c
endif
#
#
# Definition of search path variable:
#
VPATH =$(MRILU_LIB):$(MRILU_MOD):$(MRILU_OBJ)
#
# Delete target in case of errors during construction of target
#
.DELETE_ON_ERROR:
#
#
.SUFFIXES:                          # Delete the default suffixes
.SUFFIXES: .o .f .F .txt            # Define the used suffix list
#
#

#
#                    Pattern rules:
#                    ==============
#
#
include $(SOURCEDIR)/mrilu/scripts/rules.mk
#

#                        Rules implementing Actions:
#                        ===========================
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
.PHONY: install_prg
#
install_prg: default_install_prg
#
#
.PHONY: install_others
install_others:
#	cd ../getdoc;     $(MAKE) install
#	cd ../vsm;        $(MAKE) install
#	cd ../mtlb;       $(MAKE) install
#	cd ../iosrclib;   $(MAKE) install
	cd ../convprogs;  $(MAKE) install
#	cd ../misc;       $(MAKE) install
#	cd ../mtstor;     $(MAKE) install
	cd ../precon;     $(MAKE) install
#
#
.PHONY: uninstall
#
uninstall: default_uninstall

.PHONY: clean
#
clean: default_clean
	$(RM) data/Breord?? data/SchurComp?? data/ena?b.mtr data/ena?b.rhs
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
#                        Other Rules:
#                        ============
#
#
# Explicit dependencies for included files:
#
m_cmpsolprc.mod:	m_build.mod	m_cmpprc.mod	m_cmpsol.mod	m_prepro.mod\
			m_solprc.mod	m_wfree.mod	m_prcpars.mod	m_applprc.mod
m_arthura.mod:		m_build.mod	m_filgrd.mod	m_wacsr.mod
m_eugen2d.mod:      	m_build.mod	m_filgrd.mod	m_wacsr.mod
m_readsherman.mod:	m_build.mod	m_csr2csc.mod	m_fstrlen.mod	m_getunit.mod\
			m_ioerrmsg.mod	m_rdfilenm.mod	m_wacsc.mod\
			m_wfree.mod	m_guerrmsg.mod
cmpsolprc2.mod:m_build.mod	m_cmpprc.mod	m_cmpsol.mod	m_prepro.mod m_chmat.mod m_matvec.mod\
			m_solprc.mod	m_wfree.mod	m_prcpars.mod	m_applprc.mod m_rdmt.mod
#
nonsym.o:       	m_build.mod	m_cmpsolprc.mod	m_eugen2d.mod	m_fstrlen.mod	m_getunit.mod\
			m_guerrmsg.mod	m_iniglb.mod	m_matvec.mod	m_rdfilenm.mod	m_readpars.mod\
			m_rdmt.mod	m_wrdstc.mod	m_wrdvec.mod    $(mkl95_blas)
smallnonsym.o: 		m_matvec.mod	m_eugen2d.mod	m_wacsr.mod	m_iniglb.mod 	$(mkl95_blas)
nschmat.o: cmpsolprc2.mod	
#
# End of Makefile
