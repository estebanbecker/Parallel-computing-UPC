#include<upc_relaxed.h> 
#define N 100 

shared int v1[N], v2[N], v1plusv2[N]; 
void main() 
{ 
    int i; 
    for(i=MYTHREAD; i<N; i+=THREADS){
        v1plusv2[i]=v1[i]+v2[i] ; 
    }
}