# declare the compiler
CXX = g++

# declare the compiler flags flags
# -O3 flag  - optimization more for code size and execution time
CXXFLAGS  = -O3 -DAE_CPU=AE_INTEL -fopenmp -I${CURDIR}/alglib/src -std=c++11 

# other flags
LFLAGS = -L${CURDIR}/alglib/src

# libraries
LIBS=-lalglib

all: calculation

calculation: calculation.o functions_init.o
	$(CXX) $(CXXFLAGS) calculation.o functions_init.o -o  calculation $(LFLAGS) $(LIBS)

calculation.o: calculation.cpp
	$(CXX) $(CXXFLAGS) $(LFLAGS) $(LIBS) -c -Wall calculation.cpp

functions_init.o: functions_init.cpp
	$(CXX) $(CXXFLAGS) $(LFLAGS) $(LIBS) -c -Wall functions_init.cpp


clean:
	rm -rf *o calculation
