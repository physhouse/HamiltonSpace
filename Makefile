OPT     = -std=c++0x -O2
CC      = mpiicpc $(OPT)

hs_mpi: main.o Atom.o Communicator.o HamiltonSpace.o RCBTree.o
	$(CC) -o hs_mpi main.o Atom.o Communicator.o HamiltonSpace.o RCBTree.o

main.o: main.cpp
	$(CC) -c main.cpp

HamiltonSpace.o: HamiltonSpace.cpp
	$(CC) -c HamiltonSpace.cpp

Atom.o: Atom.cpp
	$(CC) -c Atom.cpp

Communicator.o: Communicator.cpp
	$(CC) -c Communicator.cpp

RCBTree.o: RCBTree.cpp
	$(CC) -c RCBTree.cpp

