CC = gcc
CFLAGS = -g -O2 -Wall

purityEM: purityEM.o EM.o read.o gamma.o
	$(CC) $(CFLAGS) purityEM.o EM.o gamma.o read.o -o purityEM -lm

purityEM.o: EM.h read.h gamma.h
	$(CC) $(CFLAGS) -c purityEM.c EM.c read.c gamma.c -lm

EM.o: EM.h
	$(CC) $(CFLAGS) -c EM.c -lm
gamma.o: gamma.h
	$(CC) $(CFLAGS) -c gamma.c -lm

read.o: read.h
	$(CC) $(CFLAGS) -c read.c -lm

clean:
	rm -rf *.o

