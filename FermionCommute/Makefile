all: build/FermionCommute ../commutators XP std

build/FermionCommute: build
	cmake -B build -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ .
	cmake --build build

XP: build/FermionCommute ../commutators
	./build/FermionCommute XP hubbard
	./build/FermionCommute XP hubbard_dispersions
	./build/FermionCommute XP continuum

std: build/FermionCommute ../commutators
	./build/FermionCommute std hubbard
	./build/FermionCommute std hubbard_dispersions
	./build/FermionCommute std continuum

continuum_xp: build/FermionCommute ../commutators
	./build/FermionCommute XP continuum

test: build/FermionCommute
	./build/FermionCommute test

debug: build/FermionCommute
	./build/FermionCommute debug

build:
	mkdir -p build

../commutators:
	mkdir -p ../commutators

clean:
	rm -rf build
	rm -rf ../commutators/*
