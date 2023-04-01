#
#
# -DUSE_DOUBLE_PRECISION for double precision (the most efficient option)
# -DSSE2 -O3 for pentium4 optimizations
# -DPARALLEL is you have mpi then compile with mpiCC and run with mpirun 
#
# NOTE: I hate makefiles and I do not use them. 
# If you know how to write a better makefile send it to me please!
#

CXX     = g++ -I../Libraries 
# Basic compilation flags (-O2 for speed)
CXXFLAGS = -O3
# Flags for double precision and sse (P4 optiomizations)
#CXXFLAGS = -DSSE2 -O3
# Flags for parallel (mpi) double precision and sse
#CXXFLAGS = -DPARALLEL -DUSE_DOUBLE_PRECISION -DSSE2 -O2

all:
	${CXX} ${CXXFLAGS} example01.cpp -o example01.exe
	${CXX} ${CXXFLAGS} example02.cpp -o example02.exe
	${CXX} ${CXXFLAGS} example03.cpp -o example03.exe
	${CXX} ${CXXFLAGS} example04.cpp -o example04.exe
	${CXX} ${CXXFLAGS} example05.cpp -o example05.exe
	${CXX} ${CXXFLAGS} example06.cpp -o example06.exe
	${CXX} ${CXXFLAGS} example07.cpp -o example07.exe
	${CXX} ${CXXFLAGS} example08.cpp -o example08.exe
	${CXX} ${CXXFLAGS} example09.cpp -o example09.exe
	${CXX} ${CXXFLAGS} example10.cpp -o example10.exe
	${CXX} ${CXXFLAGS} example11.cpp -o example11.exe
	${CXX} ${CXXFLAGS} example12.cpp -o example12.exe
	${CXX} ${CXXFLAGS} example13.cpp -o example13.exe
	${CXX} ${CXXFLAGS} example14.cpp -o example14.exe
	${CXX} ${CXXFLAGS} example15.cpp -o example15.exe
	${CXX} ${CXXFLAGS} test_sse2.cpp -o test_sse2.exe
	${CXX} ${CXXFLAGS} test_actions.cpp -o test_actions.exe
	${CXX} ${CXXFLAGS} gauge.cpp        -o gauge.exe
	${CXX} ${CXXFLAGS} meson.cpp        -o meson.exe
	${CXX} ${CXXFLAGS} ks_meson.cpp     -o ks_meson.exe
	${CXX} ${CXXFLAGS} make_meson.cpp   -o make_meson.exe
	${CXX} ${CXXFLAGS} make_meson2.cpp  -o make_meson2.exe

test:
	./test_actions.exe
	./gauge.exe -hot -nc=3 -L=12x4x4x4 -output=gauge00.mdp
	./gauge.exe -input=gauge00.mdp -beta=6.0 -steps=10 -output=gauge01.mdp
	./meson.exe -input=gauge01.mdp -lkappa=0.1 -hkappa=0.1 -ps -clean
	./ks_meson.exe -input=gauge01.mdp -lmass=0.1 -hmass=0.1 -u0=1.0 -meson=5x5


