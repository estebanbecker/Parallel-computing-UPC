#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <upc_relaxed.h>
#include <upc_collective.h> 

#define N 6

shared [(N+2)*(N+2)/THREADS] double grid[N+2][N+2], new_grid[N+2][N+2];
shared [(N+2)/THREADS] double *shared ptr[N+2], *shared new_ptr[N+2];
shared double dTmax[THREADS];
shared double diffmax;
shared double *ptr_priv[(N+2)/THREADS];
shared double *new_ptr_priv[(N+2)/THREADS];

void initialize(void)
{
    int j;

    /* Heat one side of the solid */
    for(j=1;j<N+1; j++)
    {
        grid[0][j] = 1.0;
        new_grid[0][j] = 1.0;
    }
    for( j=MYTHREAD*((N+2)/THREADS); j<(MYTHREAD+1)*((N+2)/THREADS); j++)
    {
        ptr[j] = &grid[j][0];
        new_ptr[j] = &new_grid[j][0];
    }
    int i=0;
    upc_barrier;
    for( j=MYTHREAD*((N+2)/THREADS); j<(MYTHREAD+1)*((N+2)/THREADS); j++)
    {
        ptr_priv[i] = ptr[j];
        new_ptr_priv[i] = new_ptr[j];
        i++;
    }

    upc_barrier;
}

int main(void)
{
    struct timeval ts_st, ts_end;
    double dT, epsilon, time;
    int finished, i, j, k, l;
    double T;
    int nr_iter;

    initialize();
    printf("Initialisation done\n");
    print_grid();
    /* Set the precision wanted */
    epsilon  = 0.0001;
    finished = 0;
    nr_iter = 0;

    /* and start the timed section */
    gettimeofday( &ts_st, NULL );

    do
    {
        dTmax[MYTHREAD] = 0.0;
        if(MYTHREAD!=0)
        {
            for( j=1; j<N+1; j++ )
            {  
                int i=0;
                T = 0.25 *
                    (ptr_priv[i+1][j] + ptr[(MYTHREAD*((N+2)/THREADS))-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                dT = T - ptr_priv[i][j]; /* local variation */
                if( dT > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = dT;
                new_ptr_priv[i][j] = T;
            }
        }
        
        upc_barrier;
        printf("Barrier 1 done\n");
        print_grid();

        for( i=1; i<(N+2)/THREADS-1; i++)
        {
            for( j=1; j<N+1; j++ )
            {
                T = 0.25 *
                    (ptr_priv[i+1][j] + ptr_priv[i-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                dT = T - ptr_priv[i][j]; /* local variation */
                if( dT > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = dT;
                new_ptr_priv[i][j] = T;
            }
        }
        upc_barrier;
        printf("Barrier 2 done\n");
        print_grid();

        if(MYTHREAD!=THREADS-1)
        {
            for( j=1; j<N+1; j++ )
            {
                i=(N+2)/THREADS-1;
                T = 0.25 *
                    (ptr[N+2][j] + ptr_priv[i-1][j] +
                     ptr_priv[i][j-1] + ptr_priv[i][j+1]); /* stencil */
                dT = T - ptr_priv[i][j]; /* local variation */
                if( dT > dTmax[MYTHREAD] )
                    dTmax[MYTHREAD] = dT;
                new_ptr_priv[i][j] = T;
            }
        }
        upc_barrier;
        upc_all_reduceD( &diffmax, dTmax, UPC_MAX,THREADS,1,NULL,UPC_OUT_ALLSYNC);
        if( diffmax < epsilon ) /* is the precision reached good enough ? */
            finished = 1;
        else
        {
            double shared *tmp;
            shared double *tmp_priv;
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
        if(MYTHREAD==0)
        {
            printf("Iteration %d, diffmax = %f\n", nr_iter, diffmax);
        }
        printf("Barrier 3 done\n");
        print_grid();
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
    if(MYTHREAD==0)
    {
        printf("Grid:\n");
        int i, j;
        for( i=0; i<N+2; i++ )
        {
            for( j=0; j<N+2; j++ )
                printf("%f ", new_grid[i][j]);
            printf("\n");
        }
        scanf("%d", &i);
    }
    
    upc_barrier;
    return 0;

}