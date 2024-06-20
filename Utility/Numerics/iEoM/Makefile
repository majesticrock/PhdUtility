CXX = g++

# -isystem rather than -I as this supresses warnings that occur during
# the compilation of the eigen library (I cant fix them anyways)
INCLUDEFLAGS = -isystem ~/usr/local/include

CXXFLAGS = -g $(WARNINGS) -std=c++17 $(OPT) -fopenmp

LDLIBS = -L/sw/gcc/5.3.0/rtf/lib64 -L/home/joshua/usr/local/include/boost_lib/ -lboost_serialization

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3

SRCS=ieom_test.cpp

OBJS=$(addprefix build/, $(subst .cpp,.o,$(SRCS)))

all: build build/main
	./build/main

build/main: $(OBJS) | build
	$(CXX) $(INCLUDEFLAGS) -o build/main $(OBJS) $(CXXFLAGS) $(LDLIBS)

build/%.o: %.cpp
	$(CXX) $(INCLUDEFLAGS) $< -o $@ -c $(CXXFLAGS)

build:
	mkdir -p build

clean:
	rm -rf build