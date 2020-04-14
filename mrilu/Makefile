#                                                            2002-10-10
#
#                        Description:
#                        ============
#
# Makefile for the complete installation of the programs, subroutines
# and documentation for the MRILU preconditioner software.
#
# The environment variables  NUMBASE  and  OSTYPE  should be defined,
# either with the command:
#    . ./usenumath   # to use the files in the subdirectories of the  user
#                    # and group  numath. If the modified MRILU software
#                    # should be incorporated into the CVS database then
#                    # this command should be preceded with the command:
#    . ./chgnumath   # to change the group ID into  numath
# or the command:
#    . ./lusenumath  # to update/use the files in the directory ~/numath.
#
# This make file assumes that the following subdirectories have been
# extracted from the CVS database:
#    convprogs, getdoc, iosrclib, misc, mtlb, precon, vsm
#
# This make file assumes the existence of the following subdirectories
# of the directory  $(NUMBASE):
#    bin, dat, doc, include, lib, mtlb
# These files can be created using the command:
#    ./create_numath_dirs
#
# This make file assumes the existence of the following make include
# files in the directory  $(NUMBASE)/include:
#    c_comp.mk  and  fortran_comp.mk
# The command
#    ./copy_make_inc_files
# copies the local '.mk' files, c_comp.mk  and  fortran_comp.mk, to the
# directory  $(NUMBASE)/include.
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
#   environment variable  OSTYPE  and/or the value of the options
#   CC or F77.
#   The preprocessor- and compilation options depend on the value of the
#   OPTion  variable definition, OPT, and the compilers to be used.
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
#     The name of C compiler is <name C-compiler>.  If this option is
#     omitted the name of the C compiler depends on the value of the
#     OSTYPE  environment variable.
#   F77=<name Fortran 77-compiler>
#     The name of Fortran 77 compiler is <name Fortran 77-compiler>.  If this
#     option is omitted the name of the Fortran 77 compiler depends on the
#     value of the  OSTYPE  environment variable.
#
#   The following actions are defined:
#   clean  (the default action)
#     remove all redundant files from the subdirectories.
#   install
#   . place the executable files of the changed source files in the directory
#     $(NUMBASE)/bin/
#   . place the document files of the changed source files in the directory
#     $(NUMBASE)/doc/
#   . place the object files of the changed source files in the libraries
#     $(NUMBASE)/lib/libgen.a  and  $(NUMBASE)/lib/libprecon.a
#   . place the MATLAB m-files of the changed source files in the directory
#     $(NUMBASE)/mtlb
#
#
# N.B.
# The make file uses the following environment variables:
#   NUMBASE  The base of the subdirectories in which binaries, libraries
#            etc. reside.
#   OSTYPE   The name of the operating system, which should be one of:
#            HP-UX, IRIX, IRIX64, linux, OSF1(Tru64)
#
#
#
# Define the shell 
#
SHELL = /bin/bash
#
# Delete target in case of errors during construction of target
#
.DELETE_ON_ERROR:
#
#
.SUFFIXES:                          # Delete the default suffixes
.SUFFIXES: .o .f .F .txt            # Define the used suffix list
#
# Definition of machine and user dependent variables
include makefile.inc

# Definition of search path variable:
#
VPATH = $(NUMBASE)/mrilu/scripts
##include $(NUMBASE)/mrilu/scripts/rules.mk
#
# Definition of the system dependent variables:
#   MAKE
#
# The default value:
MAKE=make
#
ifeq ($(OSTYPE), IRIX)
MAKE=gmake
endif  # IRIX
ifeq ($(OSTYPE), IRIX64)
MAKE=gmake
endif  # IRIX64
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
#
#
#                        Rules implementing Actions:
#                        ===========================
#
ifdef LIB_X11
directories=getdoc mtstor misc iosrclib convprogs precon testMRILU vsm
else
directories=getdoc mtstor misc iosrclib convprogs precon testMRILU
endif
#
# Default goal:
#
.PHONY: clean
#
clean:
	for dir in $(directories); \
	do \
	    cd $$dir; $(MAKE) clean; cd ..; \
	done
#
#
#
.PHONY: install
#
install:
	@echo Building $(directories)
	for dir in $(directories); \
	do \
	    cd $$dir; $(MAKE) install; cd ..; \
	done
#
#
.PHONY: uninstall
#
uninstall:
	for dir in $(directories); \
	do \
	    cd $$dir; $(MAKE) uninstall; cd ..; \
	done
#
#
.PHONY: test
#
test:
	cd testMRILU;	./testall $(MRILU_BIN)

# End of Makefile
