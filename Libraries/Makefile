# ifndef CXX
# CXX     = g++
# endif
# Basic compilation flags (-O2 for speed)
# Flags for double precision and sse (P4 optiomizations)
#CXXFLAGS = -DSSE2 -O3
# Flags for parallel (mpi) double precision and sse
#CXXFLAGS = -DPARALLEL -DUSE_DOUBLE_PRECISION -DSSE2 -O2

LDFLAGS := 
SOURCES_CXX := libfermiqcd.cpp
MDP_SOURCES_CXX := libmdp.cpp
CXXFLAGS_INC := -I./Libraries
CXXFLAGS := -g -DUSE_DOUBLE_PRECISION -O3 -std=c++17 -Wall -Werror -Wno-error=switch -Wno-error=sign-compare -Wno-error=attributes -pedantic $(CXXFLAGS_INC)
OBJS = $(SOURCES_CXX:.cpp=.o)
OBJS_MDP = $(MDP_SOURCES_CXX:.cpp=.o)

all: libmdp.a libfermiqcd.a

libmdp.a: $(OBJS_MDP)
	ar rc $@ $(OBJS_MDP)
	ranlib $@

libfermiqcd.a: $(OBJS)
	ar rc $@ $(OBJS)
	ranlib $@

clean:
	rm -f *.o *.a
