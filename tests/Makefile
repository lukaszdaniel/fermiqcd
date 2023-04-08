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
# -DUSE_DOUBLE_PRECISION for double precision (the most efficient option)
# -DSSE2 -O3 for pentium4 optimizations
# -DPARALLEL is you have mpi then compile with mpiCC and run with mpirun 
#
# NOTE: I hate makefiles and I do not use them. 
# If you know how to write a better makefile send it to me please!
#

# Basic compilation flags (-O2 for speed)
# Flags for double precision and sse (P4 optiomizations)
#CXXFLAGS = -DSSE2 -O3
# Flags for parallel (mpi) double precision and sse
#CXXFLAGS = -DPARALLEL -DUSE_DOUBLE_PRECISION -DSSE2 -O2

LDFLAGS := -L../Libraries -lm #-lmdpFermi
CXXFLAGS_INC := -I../Libraries
CXXFLAGS := -g -DUSE_DOUBLE_PRECISION -O3 -std=c++17 -Wall -Werror -Wno-error=switch -Wno-error=sign-compare -Wno-error=attributes -pedantic $(CXXFLAGS_INC)
#
## ==================- NOTHING TO CHANGE BELOW THIS LINE ===================
##
SRCS := $(wildcard *.cpp)
PRGS := $(patsubst %.cpp,%,$(SRCS))
PRG_SUFFIX=.exe
BINS := $(patsubst %,%$(PRG_SUFFIX),$(PRGS))
## OBJS are automagically compiled by make.
OBJS := $(patsubst %,%.o,$(PRGS))
##
all : $(BINS)
##
## For clarity sake we make use of:
.SECONDEXPANSION:
OBJ = $(patsubst %$(PRG_SUFFIX),%.o,$@)
ifeq ($(PRG_SUFFIX_FLAG),0)
	BIN = $(patsubst %$(PRG_SUFFIX),%,$@)
else
	BIN = $@
endif
## Compile the executables
%$(PRG_SUFFIX) : $(OBJS)
	$(CXX) $(OBJ)  $(LDFLAGS) -o $(BIN)
##
## $(OBJS) should be automagically removed right after linking.
##
clean:
ifeq ($(PRG_SUFFIX_FLAG),0)
	$(RM) $(PRGS)
else
	$(RM) $(BINS)
endif
##
rebuild: clean all
##
## eof Generic_Multi_Main_PWD.makefile


check:
	for file in `ls *.exe` ; do \
       echo "Testing" $${file:-3} "..."; \
	   ./$${file} > $${file%.*}".out"; \
	   diff -I "^PSIM LOG" -I "^\* Process" -I "^PROCESS " -I "... Loading time" -I "... Saving time" -I "<time>" $${file%.*}".out" $${file%.*}".save"; \
	   echo "... Done"; \
	   echo; \
	done
