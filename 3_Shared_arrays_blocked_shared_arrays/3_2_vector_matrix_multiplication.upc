// vect_mat_mult.c
#include <upc.h>
shared [THREADS] int a[THREADS][THREADS] ; //Matrix A
shared int b[THREADS], c[THREADS] ; //Vector B and C
void main (void) 

{
    int i, j; 
    //Using upc_forall to share the work
    upc_forall( i = 0 ; i < THREADS ; i++; i){
    c[i] = 0;
        for ( j= 0 ; j < THREADS ; j++){
            c[i] += a[i][j]*b[j];
        }
    }
}
