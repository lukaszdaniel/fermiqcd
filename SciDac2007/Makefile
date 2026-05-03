############################################################################
# 'A Generic Makefile for Building Multiple main() Targets in $PWD'
# Author:  Robert A. Nader (2012)
# Email: naderra at some g
# Web: xiberix
############################################################################
#  The purpose of this makefile is to compile to executable all C++ source
#  files in CWD, where each .cpp file has a main() function, and each object
#  links with a common LDFLAG.
#
#  This makefile should suffice for simple projects that require building
#  similar executable targets.  For example, if your CWD build requires
#  exclusively this pattern:
#
#  cc -c $(CXXFLAGS) main_01.cpp
#  cc main_01.o $(LDFLAGS) -o main_01
#
#  cc -c $(CXXFLAGS) main_2.cpp
#  cc main_02.o $(LDFLAGS) -o main_02
#
#  etc, ... a common case when compiling the programs of some chapter,
#  then you may be interested in using this makefile.
#
#  What YOU do:
#
#  Set PRG_SUFFIX_FLAG below to either 0 or 1 to enable or disable
#  the generation of a .exe suffix on executables
#
#  Set CXXFLAGS and LDFLAGS according to your needs.
#
#  What this makefile does automagically:
#
#  Sets SRC to a list of *.cpp files in PWD using wildcard.
#  Sets PRGS BINS and OBJS using pattern substitution.
#  Compiles each individual .cpp to .o object file.
#  Links each individual .o to its corresponding executable.
#
###########################################################################
#
PRG_SUFFIX_FLAG := 1
#
#
#
# -DUSE_SINGLE_PRECISION for single precision
# -DPARALLEL is you have mpi then compile with mpiCC and run with mpirun 
#
# NOTE: I hate makefiles and I do not use them. 
# If you know how to write a better makefile send it to me please!
#

# Basic compilation flags (-O2 for speed)
# Flags for double precision
#CXXFLAGS = -O3
# Flags for parallel (mpi) single precision
#CXXFLAGS = -DPARALLEL -DUSE_SINGLE_PRECISION -O2

# ================== Configuration ==================
#LIB_DIR := $(CURDIR)/Libraries

# ================== Common rules & flags ===========
CUSTOM_CLEAN  := 1
include ../Makefile.common

# ================== Subdirectories =================
SUBDIRS := meson meson_wave qqbar2 run1 staggered_quark wilson_quark

# ================== Default target =================
all: $(BINS) $(SUBDIRS)

.PHONY: all clean rebuild plaquette topological $(SUBDIRS)

# ================== Recursive make =================
$(SUBDIRS):
	$(MAKE) -C $@

# ================== Build rules ====================
.SECONDEXPANSION:

ifeq ($(PRG_SUFFIX_FLAG),0)
    BIN = $(patsubst %$(PRG_SUFFIX),%,$@)
else
    BIN = $@
endif

# One cpp -> one binary
%$(PRG_SUFFIX): %.o
	$(CXX) $< $(LDFLAGS) -o $(BIN)

# ================== Cleaning =======================
clean:
	@for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done
ifeq ($(PRG_SUFFIX_FLAG),0)
	$(RM) $(PRGS) *.o
else
	$(RM) $(BINS) *.o
endif

rebuild: clean all


plaquette:
	python iqcd.py code +cold +loop n=100 { +heatbath +plaquette }
	python iqcd.py compile
	./latest.exe > latest.log	
	python ibootstrap.py latest.log '"plaquette"'
	python iplot.py -p quartz

topological:
	python iqcd.py code +hot +heatbath steps=100 +loop n=100 { +ape_smear steps=1 +topological_charge }
	python iqcd.py compile
	./latest.exe > latest.log
	python ibootstrap.py latest.log '"topological_charge"'
	python isee.py topological_charge_00000.vtk

