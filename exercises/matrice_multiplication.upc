// vect_mat_mult.c
#include <upc.h>
shared [THREADS] int a[THREADS][THREADS] ;
shared int b[THREADS], c[THREADS] ;
void main (void) 

{
    int i, j; 
    upc_forall( i = 0 ; i < THREADS ; i++; i){
    c[i] = 0;
        for ( j= 0 ; j < THREADS ; j++){
            c[i] += a[i][j]*b[j];
        }
    }
}
