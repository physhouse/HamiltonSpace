OPT     = -std=c++0x -O3 -g
CC      = mpiicpc $(OPT)

hs_mpi: main.o Atom.o Communicator.o HamiltonSpace.o
	$(CC) -o hs_mpi main.o Atom.o Communicator.o HamiltonSpace.o

main.o: main.cpp
	$(CC) -c main.cpp

HamiltonSpace.o: HamiltonSpace.cpp
	$(CC) -c HamiltonSpace.cpp

Atom.o: Atom.cpp
	$(CC) -c Atom.cpp

Communicator.o: Communicator.cpp
	$(CC) -c Communicator.cpp

