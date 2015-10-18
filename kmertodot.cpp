#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>


using namespace std;

int KmerToDot(char * source, char * result){

    FILE * fp = fopen(source,"r+");
    FILE * fp1 = fopen(result,"w+");
    char * temp = new char[150000];
    while(fgets(temp, 150000, fp)){
        long int len = strlen(temp);
        for(long int i = 0; i<len;i++){
            if(temp[i] == ' '){
                temp[i]='\0';
                break;
            }
            temp[i] = tolower(temp[i]);
        }
        fprintf(fp1,"%s;\n",temp); 
        fflush(fp1);
    }
    fclose(fp1);
}

int main(int argc, char *argv[])
{   
    KmerToDot(argv[1],argv[2]);
}
