#                                                            2002-10-07
#
# Makefile for the conversion programs.
# Make can be used to:
# . produce or update the executable files in  $(MRILU_BIN)/
# . produce or update the documentation files in  $(NUMBASE)/doc/
#
#
# Usage:
#
#   make [ OPT = [[ debug | optim | prof ] ]
#        [ CC = <name C-Compiler> ] [ F77 = <name Fortran 77-Compiler> ]
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
#     files.  The range checking and profiling is disabled.
#   OPT=debug
#     The preprocessor symbol DEBUG is defined.
#     The compiler optimisation option is disabled and the options for
#     range checking and profiling are enabled to produce the object
#     files.  Extra trace code is inserted.
#   OPT=optim
#     The preprocessor symbol OPTIM is defined.
#     The compiler optimisation options are enabled to produce the object
#     files.  Range checking and profiling are disabled.
#   OPT=prof
#     The profiling is enabled.  The compiler optimisation option is enabled,
#     if possible, and the range checking is disabled.
#     The internal check of the consistency of the administration data in
#     the WS buffer is complete.
#   CC=<name C-compiler>
#     The name of C compiler is <name C-compiler>.
#   F77=<name Fortran 77-compiler>
#     The name of Fortran 77 compiler is <name Fortran 77-compiler>.
#
#   The following actions are defined:
#   all  (the default action)
#     produce the executable programs.
#   clean
#     remove all redundant files.
#   clean_doc
#     remove all the document files (extension .txt) from the current
#     directory.
#   clean_obj
#     remove all the object files from the current directory.
#   doc
#     extract the documentation from the source files.
#   install
#     update the document files in       $(NUMBASE)/doc/,
#            the executable programs in  $(MRILU_BIN)/
#   install_doc
#     place the document files of the changed source files in the directory
#     $(NUMBASE)/doc/
#   install_prg
#     place the executable programs of the changed source files in the
#     directory $(MRILU_BIN)
#
#
# N.B.
# The make file uses the following environment variables:
#   NUMBASE  The base of the subdirectories in which binaries, libraries
#            etc. reside.

#
#                   Definitions:
#                   ============
#
#
# include userdefined constants
include ../makefile.inc
#
#
# Define the basenames of the source files of the programs:
#
progs = amt2bmt  avc2bvc  bmt2amt  bvc2avc  stc2acsr stc2bcsr vsm2bmt vsmb2bmt avczn2avc hbf2amt
#
#
#
# Load from libraries:
#
LD       = $(F90)
LDLIBDIR = $(MRILU_LIB)
LDLIBES  = iosrc mtstor misc 
#
#
#
# Definition of search path variable:
VPATH = $(MRILU_LIB):$(MRILU_MOD):$(MRILU_OBJ)
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

#                    Pattern rules:
#                    ==============
#
#
include $(SOURCEDIR)/mrilu/scripts/rules.mk
#
#
# current directory:
#
#%.txt: %.F
#	$(MRILU_BIN)/getdoc $< $@
#
#%.txt: %.f
#	$(MRILU_BIN)/getdoc $< $@
#
#
#
# Pattern rule to extract the documentation file and to place it in the
# directory  $(NUMBASE)/doc/:
#
#$(NUMBASE)/doc/%.txt: %.F
#	$(RM) $@
#	$(MRILU_BIN)/getdoc $< $@
#	-chmod u=rw,g=rw,o=r $@
#
#$(NUMBASE)/doc/%.txt: %.f
#	$(RM) $@
#	$(MRILU_BIN)/getdoc $< $@
#	-chmod u=rw,g=rw,o=r $@
#
#
#
# Default goal:
#
#.PHONY: all
#all: $(addprefix $(MRILU_BIN)/, $(progs))
#
#
#
#                        Rules implementing Actions:
#                        ===========================
#
# Default goal: compile all modules
#
#.PHONY: default_all
#
#default_all: $(addsuffix .mod, $(addprefix m_, $(modules) ) ) $(addsuffix .o, $(progs) )
#
# Installation
#
#.PHONY: default_install
#
#default_install: install_others install_doc all install_obj install_prg
#
# Making the documentation
#
#.PHONY: default_install_doc
#
#default_install_doc:\
#             $(addprefix $(NUMBASE)/doc/, $(addsuffix .txt, $(modules) ) ) \
#             $(addprefix $(NUMBASE)/doc/, $(addsuffix .txt, $(progs) ) )
#
#
# Adding the object files to the library
#
#.PHONY: default_install_obj
#
#default_install_obj:\
#             ($MRILU_LIB)/$(thislib).a($(addsuffix .o, $(modules) )) \
#
#.PHONY: default_install_prg
#
#install_prg:\
#             $(addprefix $(MRILU_BIN)/, $(progs) )\
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
doc: $(addsuffix .txt, $(progs) )
#
#
#
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
install_obj:
#
.PHONY: install_prg
#
install_prg: default_install_prg
#
#
.PHONY: install_others
install_others:
#	cd ../getdoc;     $(MAKE) install
	cd ../mtstor;     $(MAKE) install
	cd ../iosrclib;   $(MAKE) install
#
#
#
#
#
#
#                        Other Rules:
#                        ============
#

amt2bmt.o: 	m_build.mod	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtmt.mod  	m_rdmt.mod\
		m_fstrlen.mod	m_wacsr.mod	m_dump.mod
hbf2amt.o: m_build.mod	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtmt.mod  	m_rdhbf.mod\
		m_fstrlen.mod	m_wacsr.mod
bmt2amt.o: 	m_build.mod	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtmt.mod  	m_rdmt.mod\
		m_fstrlen.mod	m_wacsr.mod
avc2bvc.o: 	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtvec.mod 	m_rdvec.mod 	m_fstrlen.mod
bvc2avc.o: 	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtvec.mod 	m_rdvec.mod	 m_fstrlen.mod
avczn2avc.o: m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtvec.mod  m_rdvc.mod 	m_fstrlen.mod
stc2acsr.o:	m_build.mod	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtvec.mod 	m_wrtmt.mod\
		m_rdstc.mod	m_fstrlen.mod  	m_getunit.mod	m_wacsr.mod
stc2bcsr.o:	m_build.mod	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtvec.mod 	m_wrtmt.mod\
	 	m_rdstc.mod	m_fstrlen.mod	m_wacsr.mod
vsm2bmt.o: 	m_build.mod	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtmt.mod  	m_rdvsm.mod\
		m_fstrlen.mod	m_wacsr.mod	m_dump.mod
vsmb2bmt.o: m_build.mod	m_rdfilenm.mod 	m_ioerrmsg.mod 	m_wrtmt.mod m_rdvsmb.mod\
		m_fstrlen.mod	m_wacsr.mod	m_dump.mod
#
# End of Makefile

















