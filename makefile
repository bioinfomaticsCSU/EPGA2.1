CC=g++

CPPFLAGS = -g -Wall -O3

LDFLAGS=-lpthread

EPGA:	main.o
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp bitarray.h common.h constructcontigset.h dfs.h graph.h kmerset.h readset.h contigmerge.h fillgap.h scaffolding.h
	$(CC) -o $@ -c $< $(LDFLAGS)

GetKmerHash:	getkmerhashtable.o
	$(CC) -o $@ $^ $(LDFLAGS)

getkmerhashtable.o: getkmerhashtable.cpp bitarray.h common.h kmerset.h readset.h
	$(CC) -o $@ -c $< $(LDFLAGS)

SimplePathToGraph:	simplepathtograph.o
	$(CC) -o $@ $^ $(LDFLAGS)

simplepathtograph.o: simplepathtograph.cpp bitarray.h common.h kmerset.h readsetSimplePath.h graphSimplePath.h
	$(CC) -o $@ -c $< $(LDFLAGS)

KmerToDot:	kmertodot.o
	$(CC) -o $@ $^ $(LDFLAGS)

kmertodot.o: kmertodot.cpp
	$(CC) -o $@ -c $< $(LDFLAGS)

all: EPGA GetKmerHash SimplePathToGraph KmerToDot
	make -C Bcalm/
	cp Bcalm/bcalm ./
	make -C Bless/
	cp Bless/bless ./
	make -C Dsk/
	cp Dsk/dsk ./
	cp Dsk/parse_results ./

clean:
	rm -f *.o
	rm EPGA dsk bless bcalm parse_results GetKmerHash SimplePathToGraph KmerToDot
