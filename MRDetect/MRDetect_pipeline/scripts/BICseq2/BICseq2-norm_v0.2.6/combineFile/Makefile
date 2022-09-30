CC = gcc
CFLAGS = -g -O2 -Wall

combineFile: combine.o read.o
	$(CC) $(CFLAGS) combine.o read.o -o combineFile -lm

combine.o: read.h
	$(CC) $(CFLAGS) -c combine.c read.c -lm

read.o: read.h
	$(CC) $(CFLAGS) -c read.c -lm

clean:
	rm -rf *.o

