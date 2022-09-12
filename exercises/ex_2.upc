#include <upc_relaxed.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TOTALSIZE 800

void init();

shared double x_new[TOTALSIZE];
shared double x[TOTALSIZE];
shared double b[TOTALSIZE];

int main(int argc, char **argv){
    int j;

    init();
    upc_barrier;

    upc_forall( j=1; j<TOTALSIZE-1; j++ ; j){
        if( MYTHREAD ==  j % THREADS ){
            x_new[j] = 0.5 * ( x[j-1] + x[j+1] + b[j] );
        }
    }
    
    if( MYTHREAD == 0 ){
        printf("   b   |    x   | x_new\n");
        printf("=============================\n");
    }

    upc_forall( j=0; j<TOTALSIZE; j++ ; j)
        if( MYTHREAD ==  j % THREADS ){
            printf("%1.4f | %1.4f | %1.4f \n", b[j], x[j], x_new[j]);
        }
            

    return 0;
}

void init(){
    int i;

    if(MYTHREAD == 0){
        srand( time(NULL) );

        for( i=0; i<TOTALSIZE; i++ ){
            b[i] = (double)rand() / RAND_MAX;
            x[i] = (double)rand() / RAND_MAX;
        }
    }
}

