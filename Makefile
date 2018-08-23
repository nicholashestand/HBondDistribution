# simple makefile for order program

all: HBondDistribution.exe

HBondDistribution.exe: HBondDistribution.cpp
	g++ HBondDistribution.cpp -o HBondDistribution.exe -lm -lgmx_reader -lxdrfile -std=c++11 -fmax-errors=10

clean:
	rm -f HBondDistribution.exe
