#include<upc_relaxed.h> 
#define N 100 

shared int v1[N], v2[N], v1plusv2[N]; 
void main() 
{ 
    int i; 
    for(i=0;i<N;i++){
    if(MYTHREAD==i%THREADS) 
        v1plusv2[i]=v1[i]+v2[i] ; 
    }
}