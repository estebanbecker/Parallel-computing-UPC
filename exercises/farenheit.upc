#include <stdio.h>
#include <upc.h>
#define TBL_SZ 12

int main(){
    static shared int fahrenheit[TBL_SZ];
    static shared int step=10; 
    int celsius, i;

    upc_forall(i=MYTHREAD;i<TBL_SZ;i+=THREADS;i){
        celsius=step*i;
        fahrenheit[i]=celsius*(9.0/5.0)+32;
    }
    upc_barrier;

    if(MYTHREAD==0)

    for(i=0;i<TBL_SZ;i++){
        celsius=step*i;
        printf("%d \t %d \n", fahrenheit[i], celsius);
    }
}
