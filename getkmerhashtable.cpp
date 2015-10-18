#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <pthread.h>
#include <time.h>
#include "bitarray.h"
#include "readset.h"
#include "kmerset.h"

using namespace std;

int main(int argc, char *argv[])
{   
    int i = 0;
    long int setNumber = argc-6;
    long int kmerLength = atoi(argv[argc-4]);
    long int threadNumber = atoi(argv[argc-3]);
    unsigned long long int kmerSetHashTableCount = atoi(argv[argc-5]);
    KmerSet * kmerSet = new KmerSet[setNumber];
    for(i = 0; i<setNumber;i++){  
        kmerSet[i].address = argv[i+1];
    }
	KmerSetHashTableHead * kmerSetHashTableHead = new KmerSetHashTableHead;
    InitKmerSet(kmerSet, kmerSetHashTableHead, kmerLength, kmerSetHashTableCount, setNumber, threadNumber, argv[argc-2], argv[argc -1]);
}
