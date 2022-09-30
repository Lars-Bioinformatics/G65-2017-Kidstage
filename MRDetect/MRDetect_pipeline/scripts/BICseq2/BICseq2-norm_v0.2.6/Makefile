CC = gcc
CFLAGS = -o -g -O2 -Wall
CFLAGSEND = -lm

all: combineFile PrepPois PrepPoisGAM BRS filterCNV purityEM random_split

combineFile: combineFile/combine.o lib/read.o
	$(CC) $(CFLAGS) combineFile/combine.o lib/read.o -o combineFile/combineFile $(CFLAGSEND)

PrepPois: DataPrepare/DtaPrep.o lib/read.o lib/statlib.o
	$(CC) $(CFLAGS)  DataPrepare/DtaPrep.o lib/read.o lib/statlib.o -o DataPrepare/PrepPois $(CFLAGSEND)

PrepPoisGAM: DataPrepareGAM/DtaPrepGAM.o lib/read.o lib/statlib.o
	$(CC) $(CFLAGS) DataPrepareGAM/DtaPrepGAM.o lib/read.o lib/statlib.o -o DataPrepareGAM/PrepPoisGAM $(CFLAGSEND)

BRS: BRS_1bp/BRS.o BRS_1bp/pos_cnt_lst.o lib/read.o
	$(CC) $(CFLAGS) BRS_1bp/BRS.o BRS_1bp/pos_cnt_lst.o lib/read.o -o BRS_1bp/BRS $(CFLAGSEND)

filterCNV: Filter/filter.o lib/read.o
	$(CC) $(CFLAGS) Filter/filter.o lib/read.o -o Filter/filterCNV $(CFLAGSEND)

purityEM: purityEM/purityEM.o purityEM/EM.o lib/read.o purityEM/gamma.o
	 $(CC) $(CFLAGS) purityEM/purityEM.o purityEM/EM.o lib/read.o purityEM/gamma.o -o  purityEM/purityEM $(CFLAGSEND)

random_split: randomSample/split.o
	$(CC) $(CFLAGS) randomSample/split.o -o randomSample/random_split $(CFLAGSEND)


combine.o: combineFile/read.h
	 $(CC) $(CFLAGS) -c combineFile/combine.c lib/read.c $(CFLAGSEND)

DtaPrep.o: lib/read.h lib/statlib.h
	 $(CC) $(CFLAGS) -c DataPrepare/DtaPrep.c lib/read.c lib/statlib.c $(CFLAGSEND)

PrepPoisGAM.o: lib/read.h lib/statlib.h
	$(CC) $(CFLAGS) -c DataPrepareGAM/DtaPrepGAM.c lib/read.c lib/statlib.c $(CFLAGSEND)

BRS.o: BRS_1bp/pos_cnt_lst.h lib/read.h
	$(CC) $(CFLAGS) -c BRS_1bp/pos_cnt_lst.c lib/read.c $(CFLAGSEND)

filter.o: lib/read.h
	$(CC) $(CFLAGS) -c Filter/filter.c lib/read.c $(CFLAGSEND)


purityEM.o: purityEM/EM.h lib/read.h purityEM/gamma.h
	$(CC) $(CFLAGS) -c purityEM/purityEM.c purityEM/EM.c lib/read.c purityEM/gamma.c $(CFLAGSEND)

split.o:
	 $(CC) $(CFLAGS) -c randomSample/split.c $(CFLAGSEND)



EM.o: purityEM/EM.h
	 $(CC) $(CFLAGS) -c purityEM/EM.c $(CFLAGSEND)

gamma.o: purityEM/gamma.h
	$(CC) $(CFLAGS) -c purityEM/gamma.c $(CFLAGSEND)


pos_cnt_lst.o: BRS_1bp/pos_cnt_lst.h
	$(CC) $(CFLAGS) -c BRS_1bp/pos_cnt_lst.c $(CFLAGSEND)


statlib.o: lib/statlib.c
	$(CC) $(CFLAGS) -c lib/statlib.c $(CFLAGSEND)

read.o: lib/read.h
	$(CC) $(CFLAGS) -c lib/read.c $(CFLAGSEND)

clean:
	rm -rf combineFile/*.o DataPrepare/*.o lib/*.o DataPrepareGAM/*.o BRS_1bp/*.o Filter/*.o purityEM/*.o randomSample/*.o

