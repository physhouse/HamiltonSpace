OPT     = -std=c++0x -O3
CC      = g++ $(OPT)

hs_mpi: main.o NMDPD.o Frame.o
	$(CC) -o NMDPD.x main.o NMDPD.o Frame.o

main.o: main.cpp
	$(CC) -c main.cpp

NMDPD.o: NMDPD.cpp
	$(CC) $(CFLAGS) -c NMDPD.cpp

Frame.o: Frame.cpp
	$(CC) -c Frame.cpp


