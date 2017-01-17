#                                                            2003-02-06
#
#                        Description:
#                        ============
#
# Makefile for the matrix storage and help subroutines in  mtstor/.
# Make can be used to:
# . produce or update the documentation files in  $(NUMBASE)/doc/
# . update the include files in  $(NUMBASE)/include/
# . produce or update the object files in  $(MRILU_LIB)/libprecon.a
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
#     produce the object files.
#   clean
#     remove all redundant files.
#   clean_doc
#     remove all the document files (extension '.txt') from the current
#     directory.
#   clean_ftn
#     remove each '.f', with an existing pre-Fortran 77 '.F' file, from the
#     current directory.
#   clean_obj
#     remove all the object files from the current directory.
#   doc
#     extract the documentation from the source files.
#   fortran_src
#     produce Fortran source files, with extension '.f', from the pre-Fortran
#     source files with extension '.F'.
#   install
#     update the document files in    $(NUMBASE)/doc/,
#            the include files in     $(NUMBASE)/include  and
#            the object files in      $(MRILU_LIB)/libprecon.a
#   install_doc
#     place the document files of the changed source files in the directory
#     $(NUMBASE)/doc/
#   install_inc
#     place the include file(s) in the directory  $(NUMBASE)/include/
#   install_obj
#     place the object files of the changed source files in the library
#     $(MRILU_LIB)/libprecon.a
#
#
# N.B.
# The make file uses the following environment variables:
#   NUMBASE  The base of the subdirectories in which binaries, libraries
#            etc. reside.
#
#                        Definitions:
#                        ============
# include userdefined constants
include ../makefile.inc
#
#
# Define the name of the directory where the Fortran include files are stored:
#
FINCLDIR = $(NUMBASE)/include
#
#
# Definition of the system dependent variables for the compilers:
# the FPP preprocessor, the Fortran 77 compiler  and
# the FPP preprocessor flags, the Fortran 77 compiler flags:
#
#include $(NUMBASE)/include/fortran_comp.mk
#
#
#
# Define the base names of the source files of the subroutines involved:
#
modules =  dump  build     cscvec  csr2csc csrvec	chkcnt    \
           mterrmsg  \
           wacsc     wacsr     wacsrd    wadia     wafm      \
           wamlp     wapffp    wapldu    wappnd    wapsfp    \
           waprc     wascbm    wascde    wcompr    \
           wennz     csrresize cscresize wfree
objects =
#
# Load from libraries:
#
LDLIBDIR = $(MRILU_LIB)
LDLIBES  = -lmisc

thislib = libmtstor
#
#
# Define the file names of the include files in the current directory:
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
#
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
##	cd ../misc;       $(MAKE) install
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
.PHONY: doc
#
doc: $(addsuffix .txt, $(modules) )
#
#
#
.PHONY: fortran_src
#
fortran_src: $(addsuffix .f, $(modules) )
#

#                        Other Rules:
#                        ============
#
#
# Explicit dependencies for included files:
#
preconlib = $(MRILU_LIB)/libprecon.a
#
m_build.mod:		m_dump.mod
m_chkcnt.mod:		m_build.mod	m_mterrmsg.mod
m_mterrmsg.mod: 	m_build.mod
m_wacsc.mod:    	m_build.mod    
m_wacsr.mod:    	m_build.mod    
m_wacsrd.mod:   	m_build.mod    	m_wadia.mod	m_wacsr.mod
m_wadia.mod:    	m_build.mod    
m_wafm.mod:     	m_build.mod    
m_wamlp.mod:    	m_build.mod    
m_wapffp.mod:   	m_build.mod    
m_wapldu.mod:   	m_build.mod    	m_wacsr.mod	m_wadia.mod	m_wafm.mod
m_wappnd.mod:   	m_build.mod
m_waprc.mod:    	m_build.mod
m_wapsfp.mod:   	m_build.mod     m_wacsr.mod	m_wadia.mod
m_wascbm.mod:   	m_build.mod	m_wacsr.mod	m_wadia.mod
m_wascde.mod:   	m_build.mod	m_wacsr.mod	m_wadia.mod
m_wcompr.mod:   	m_build.mod	m_cscresize.mod	m_csrresize.mod
m_wennz.mod:    	m_build.mod
m_wextend.mod:   	m_build.mod
m_wfree.mod:    	m_build.mod
#
# End of Makefile

