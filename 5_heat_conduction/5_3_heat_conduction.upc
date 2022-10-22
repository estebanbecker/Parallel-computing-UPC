#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <upc_relaxed.h>
#include <upc_collective.h> 

#define N 6
#define PRIV_SIZE ((N+2)/THREADS)

shared [(N+2)*(N+2)/THREADS] double grid[N+2][N+2], new_grid[N+2][N+2];
shared double dTmax[THREADS];
shared double diffmax;
double *new_ptr_priv[PRIV_SIZE], *ptr_priv[PRIV_SIZE], *tmp_priv;
shared [(N+2)*(N+2)/THREADS] double *ptr[N+2], *new_ptr[N+2], *tmp;

void initialize(void)
{
    int j;

    /* Heat one side of the solid */
    for(j=1;j<N+1; j++)
    {
        grid[0][j] = 1.0;
        new_grid[0][j] = 1.0;
    }
}

int main(void)
{
    struct timeval ts_st, ts_end;
    double dT, epsilon, time;
    int finished, i, j, k, l;
    double T;
    int nr_iter;

    if(MYTHREAD == 0)
    {
        initialize();
    }
    upc_barrier;
    /* Set the precision wanted */
    epsilon  = 0.0001;
    finished = 0;
    nr_iter = 0;

    /*Initialize the pointers*/
    for (i = 0; i < N + 2; i++)
    {
        ptr[i] = grid[i];
        new_ptr[i] = new_grid[i];
    }

    for (i = 0; i < PRIV_SIZE; i++)
    {
        ptr_priv[i] = (double *) grid[i + (MYTHREAD * PRIV_SIZE)];
        new_ptr_priv[i] = (double *) new_grid[i + (MYTHREAD * PRIV_SIZE)];
    }   

    /* and start the timed section */
    gettimeofday( &ts_st, NULL );

    do
    {
        dTmax[MYTHREAD] = 0.0;
        if(MYTHREAD!=0)
        {
            /*First for loop for the first line*/
            for( j=1; j<N+1; j++ )
            {  
                i=0;
                T = 0.25 *
                    (ptr_priv[i + 1][j] + ptr[(MYTHREAD*PRIV_SIZE)-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                
                dT = T - ptr_priv[0][j]; /* local variation */
                if( fabs(dT) > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = fabs(dT);
                new_ptr_priv[i][j] = T;
            }
        }
        
        /*Second for loop for the lines in the middle with only private pointers*/
        for( i=1; i<(N+2)/THREADS-1; i++)
        {
            for( j=1; j<N+1; j++ )
            {
                T = 0.25 *
                    (ptr_priv[i+1][j] + ptr_priv[i-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                dT = T - ptr_priv[i][j]; /* local variation */
                
                if( fabs(dT) > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = fabs(dT);
    
                new_ptr_priv[i][j] = T;
            }
        }

        /*Third for loop for the last line*/
        if(MYTHREAD!=THREADS-1)
        {
            for( j=1; j<N+1; j++ )
            {
                T = 0.25 *
                    (ptr[(MYTHREAD + 1)* PRIV_SIZE][j] + ptr_priv[i-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                dT = T - ptr_priv[i][j]; /* local variation */

                if( fabs(dT) > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = fabs(dT);
                new_ptr_priv[i][j] = T;
            }
        }
        upc_all_reduceD( &diffmax, dTmax, UPC_MAX,THREADS,1,NULL,UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);

        
        if( diffmax < epsilon ) /* is the precision reached good enough ? */
            finished = 1;
        else
        {
            /*Swap the pointers*/
            upc_forall(k = 0 ; k<N+2 ; k++ ; k)   
            {   
                tmp = ptr[k];
                ptr[k] = new_ptr[k];
                new_ptr[k] = tmp;  
            }
            for( l=0; l<PRIV_SIZE; l++)
            {
                tmp_priv = ptr_priv[l];
                ptr_priv[l] = new_ptr_priv[l];
                new_ptr_priv[l] = tmp_priv;
            }
        }
        upc_barrier;
        nr_iter++;
    } while( finished == 0 );

    if(MYTHREAD == 0)
    {
        gettimeofday( &ts_end, NULL ); /* end the timed section */

        /* compute the execution time */
        time = ts_end.tv_sec + (ts_end.tv_usec / 1000000.0);
        time -= ts_st.tv_sec + (ts_st.tv_usec / 1000000.0);

        printf("%d iterations in %.5lf sec\n", nr_iter, time);
    }
    
    return 0;
}

//A function to print the grid
int print_grid()
{   
    upc_barrier;
    if(MYTHREAD==0)
    {
        printf("Grid:\n");
        int i, j;
        for( i=0; i<N+2; i++ )
        {
            for( j=0; j<N+2; j++ )
                printf("%f ", new_ptr[i][j]);
            printf("\n");
        }
        
    }
    upc_barrier;
    return 0;

}