CXX = g++
CXXFLAGS = -O3 -march=x86-64 -std=c++14 -I../Libraries -DLINUX
SOURCES = qcdsun.cpp
EXECUTABLE = qcdsun

$(EXECUTABLE): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(SOURCES) 

clean:
	rm -f *.o
