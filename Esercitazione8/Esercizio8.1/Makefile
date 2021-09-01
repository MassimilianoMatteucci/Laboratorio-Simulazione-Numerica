
CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o library.o
	$(CC) library.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
library.o : library.cpp library.h
	$(CC) -c library.cpp -o library.o $(CFLAGS)
clean :
	rm *.o main.exe 
