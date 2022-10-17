#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <upc_relaxed.h>
#include <upc_collective.h> 

#define N 4
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

    for(i=0; i<N+2; i+=THREADS)
    {
        ptr[i] = grid[i*(N+2)/THREADS];
        new_ptr[i] = new_grid[i*(N+2)/THREADS];
    }

    for (i = 0; i < N + 2; i++)
    {
        ptr[i] = &grid[i][0];
        new_ptr[i] = &new_grid[i][0];
    }

    for (i = 0; i < PRIV_SIZE; i++)
    {
        ptr_priv[i] = (double *)&grid[i + (MYTHREAD * PRIV_SIZE)][0];
        new_ptr_priv[i] = (double *)&new_grid[i + (MYTHREAD * PRIV_SIZE)][0];
    }   

    /* and start the timed section */
    gettimeofday( &ts_st, NULL );

    do
    {
        dTmax[MYTHREAD] = 0.0;
        if(MYTHREAD!=0)
        {
            for( j=1; j<N+1; j++ )
            {  
                i=0;
                T = 0.25 *
                    (ptr_priv[i + 1][j] + ptr[(MYTHREAD*((N+2)/THREADS))-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                dT = T - ptr_priv[0][j]; /* local variation */
                if( dT > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = fabs(dT);
                new_ptr_priv[i][j] = T;
            }
        }
        
        upc_barrier;

        for( i=1; i<(N+2)/THREADS-1; i++)
        {
            for( j=1; j<N+1; j++ )
            {
                T = 0.25 *
                    (ptr_priv[i+1][j] + ptr_priv[i-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                dT = T - ptr_priv[i][j]; /* local variation */
                if( dT > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = fabs(dT);
                new_ptr_priv[i][j] = T;
            }
        }
        upc_barrier;

        if(MYTHREAD!=THREADS-1)
        {
            for( j=1; j<N+1; j++ )
            {
                i=(N+2)/THREADS-1;
                T = 0.25 *
                    (ptr[(MYTHREAD + 1)* PRIV_SIZE][j] + ptr_priv[i-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                dT = T - ptr_priv[i][j]; /* local variation */
                if( dT > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = fabs(dT);
                new_ptr_priv[i][j] = T;
            }
        }
        upc_barrier;
        upc_all_reduceD( &diffmax, dTmax, UPC_MAX,THREADS,1,NULL,UPC_OUT_ALLSYNC);
        
        printf("nr_iter = %d, diffmax = %f", nr_iter, diffmax);
        print_grid();
        
        if( diffmax < epsilon ) /* is the precision reached good enough ? */
            finished = 1;
        else
        {
            int i=0;
            for( k=MYTHREAD*((N+2)/THREADS); k<(MYTHREAD+1)*((N+2)/THREADS); k++)     /* not yet ... Need to prepare */
            {                                    /* the next iteration */
                tmp = ptr[k];
                ptr[k] = new_ptr[k];
                new_ptr[k] = tmp;

                tmp_priv=ptr_priv[i];
                ptr_priv[i]=new_ptr_priv[i];
                new_ptr_priv[i]=tmp_priv;

                i++;
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

    print_grid();
    
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
        scanf("%d", &i);
    }
    
    upc_barrier;
    return 0;

}