CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o library.o
	$(CC) library.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
library.o : library.cpp library.h
	$(CC) -c library.cpp -o library.o $(CFLAGS)

clean :
	rm *.o Monte_Carlo_ISING_1D.exe
	rm -rf output.*
	rm -rf seed.out
	rm -rf config.final
