# ifndef CXX
# CXX     = g++
# endif
# Basic compilation flags (-O2 for speed)
# Flags for double precision
#CXXFLAGS = -O3
# Flags for parallel (mpi) single precision
#CXXFLAGS = -DPARALLEL -DUSE_SINGLE_PRECISION -O2

# ================== Configuration ==================
LIB_DIR := $(CURDIR)

# ================== Common flags ===================
LIBRARY_BUILD := 1
CUSTOM_CLEAN  := 1
include ../Makefile.common

# ================== Sources ========================
SOURCES_CXX     := libfermiqcd.cpp
MDP_SOURCES_CXX := libmdp.cpp

OBJS     := $(SOURCES_CXX:.cpp=.o)
OBJS_MDP := $(MDP_SOURCES_CXX:.cpp=.o)

# ================== Libraries ======================
all: libmdp.a libfermiqcd.a

libmdp.a: $(OBJS_MDP)
	$(AR) rc $@ $^
	$(RANLIB) $@

libfermiqcd.a: $(OBJS)
	$(AR) rc $@ $^
	$(RANLIB) $@

# ================== Cleaning =======================
clean:
	$(RM) *.o *.a
