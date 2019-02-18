#                                                           2002-10-07
# Makefile for the executables:
#   getdoc  and  showdoc
#
# Usage:
#
#   make [ OPT = [ debug | optim | prof ] ]
#        [ CC = <name C-Compiler> ]
#        <action>
#
#   The contents of the C source files depend on the value of the OPTion
#   variable definition, OPT.
#   The preprocessor and compiler depend on the value of the environment
#   variable  OSTYPE  and/or the value of the option  CC.
#   The compilation options depend on the value of the  OPTion  variable
#   definition, OPT, and the compiler to be used.
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
#
#   The following values of <action> are defined:
#   all      make the executable program  getdoc.
#   <empty>  make the executable program  getdoc.
#   clean    remove all object files, executable programs  and  the
#            files matching *~.
#   install  make the executable program  getdoc  and copy it to
#            $(MRILU_BIN)/.
#            Copy the script  showdoc  to  $(MRILU_BIN)/
#
#
# N.B.
# The make file uses the following environment variables:
#   NUMBASE  The base of the subdirectories in which binaries, libraries
#            etc. reside.
#   OSTYPE   The name of the operating system, which should be one of:
#            HP-UX, IRIX, IRIX64, linu, OSF1(Tru64)
#
###############################################################################
#
#                        Definitions:
#                        ============
# include userdefined constants
include ../makefile.inc
#
# Define the base names of the modules involved:
#
modules =
#
# Define the base names of the main programs involved:
#
progs = getdoc
#
# Define the name of the library to be made:
#
thislib = getdoc
#
# Define the loader and its environment
#
LD       = $(CC)
LDLIBDIR = 
LDLIBES  = 
#
#
# Definition of search path variable:
#
VPATH = $(MRILU_LIB):$(MRILU_MOD):$(MRILU_OBJ)
#
#
#
# Delete target in case of errors during construction of target
#
.DELETE_ON_ERROR:
#
#
.SUFFIXES:                          # Delete the default suffixes
.SUFFIXES: .o .c .h .F90 .txt       # Define the used suffix list
#
###############################################################################
#
#
#                    Pattern rules:
#                    ==============
#
include $(SOURCEDIR)/mrilu/scripts/rules.mk
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
install_doc: default_install_doc $(MRILU_BIN)/showdoc
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
.PHONY: doc
#
doc: $(addsuffix .txt, $(modules) ) 
#
#
$(MRILU_BIN)/showdoc:
	$(RM) $@
	echo NUMBASE=$(NUMBASE) > showdoc
	cat showdoc_template.sh >> showdoc
	cp showdoc $@
	-chmod u=rwx,g=rwx,o=rx $@
#
$(NUMBASE)/doc/getdoc.txt: getdoc.txt
	$(RM) $@
	cp  getdoc.txt $@
	-chmod u=rw,g=rw,o=r $@
#
#
#
#                        Other Rules:
#                        ============
#
# Explicit definition of dependencies:
#
gendefs.o: gendefs.h
getdoc.o:  gendefs.h
#
# End of Makefile
