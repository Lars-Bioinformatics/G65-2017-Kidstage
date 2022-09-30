CC = gcc
CFLAGS = -g -O2 -Wall



all: combineFile bootstrap EstLambdaFct countRead genotype


combineFile: combineFile/combine.o lib/read.o
	$(CC) $(CFLAGS) combineFile/combine.o lib/read.o -o combineFile/combineFile -lm


bootstrap:  bootstrap/bootstrap.o lib/read.o lib/gamma.o
	$(CC) $(CFLAGS) bootstrap/bootstrap.o lib/read.o lib/gamma.o -o bootstrap/bootstrapTest -lm


EstLambdaFct: EstimateLambdaFactor/EstLambdaFct.o lib/read.o
	 $(CC) $(CFLAGS) EstimateLambdaFactor/EstLambdaFct.o lib/read.o -o EstimateLambdaFactor/EstLambdaFct -lm
	
countRead: genotype/countRead.o lib/read.o genotype/ITree.o lib/rbtree.o
	$(CC) $(CFLAGS) genotype/countRead.o lib/read.o genotype/ITree.o lib/rbtree.o -o genotype/countRead -lm

genotype: genotype/genotype.o lib/read.o lib/gamma.o
	$(CC) $(CFLAGS) genotype/genotype.o lib/read.o lib/gamma.o -o genotype/genotype -lm
	
MBICseq: MBIC-seq_v0.1.5/bin.o MBIC-seq_v0.1.5/mbic-seq.o lib/rbtree.o lib/read.o
	 $(CC) $(CFLAGS) MBIC-seq_v0.1.5/bin.o MBIC-seq_v0.1.5/mbic-seq.o lib/rbtree.o lib/read.o -o MBIC-seq_v0.1.5/MBICseq -lm








combine.o: combineFile/read.h
	$(CC) $(CFLAGS) -c combineFile/combine.c lib/read.c -lm

bootstrap.o: lib/read.h
	$(CC) $(CFLAGS) -c bootstrap/bootstrap.c lib/read.c -lm

EstLambdaFct.o: lib/read.h
	$(CC) $(CFLAGS) -c EstimateLambdaFactor/EstLambdaFct.c lib/read.c -lm

countRead.o: lib/read.h
	$(CC) $(CFLAGS) -c genotype/countRead.c lib/read.c genotype/ITree.c lib/rbtree.c -lm
genotype.o: lib/read.h
	$(CC) $(CFLAGS) -c genotype/genotype.c lib/read.c -lm

bin.o: MBIC-seq_v0.1.5/bin.h lib/rbtree.h lib/read.h
	$(CC) $(CFLAGS) -c MBIC-seq_v0.1.5/bin.c lib/read.c lib/rbtree.c -lm





ITree.o: genotype/ITree.h
	$(CC) $(CFLAGS) -c genotype/ITree.c lib/rbtree.c -lm

rbtree.o: lib/rbtree.h
	$(CC) $(CFLAGS) -c lib/rbtree.c -lm

gamma.o: lib/gamma.h
	$(CC) $(CFLAGS) -c lib/gamma.c -lm

read.o: lib/read.h
	$(CC) $(CFLAGS) -c lib/read.c -lm


clean:
	rm -rf lib/*.o combineFile/*.o bootstrap/*.o EstimateLambdaFactor/*.o genotype/*.o MBIC-seq_v0.1.5/*.o
