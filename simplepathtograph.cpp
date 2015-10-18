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
#include "readsetSimplePath.h"
#include "kmerset.h"
#include "graphSimplePath.h"

using namespace std;

int main(int argc, char *argv[])
{   
    GetSimplePathFromBcalm(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));
}
