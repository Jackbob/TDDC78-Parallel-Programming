CC = g++
CPP = /usr/local/opt/llvm/bin/clang
CXXFLAGS = -Xpreprocessor -fopenmp /usr/local/opt/libomp/include
LDFLAGS = /usr/local/opt/libomp/lib/
example: laplsolv.cpp
	$(CC) $(CPP) $(CPPFLAGS) -o example laplsolv.cpp $(LDFLAGS)