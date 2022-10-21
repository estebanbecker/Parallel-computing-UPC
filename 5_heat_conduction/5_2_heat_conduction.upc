#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <upc_relaxed.h>
#include <upc_collective.h> 

#define N 30

shared [(N+2)*(N+2)/THREADS] double grid[N+2][N+2], new_grid[N+2][N+2];
shared [(N+2)*(N+2)/THREADS] double *shared ptr[N+2], *shared new_ptr[N+2];
shared double dTmax[THREADS];
shared double diffmax;

void initialize(void)
{
    int j;

    /* Heat one side of the solid */
    upc_forall( j=1; j<N+1; j++; j)
    {
        grid[0][j] = 1.0;
        new_grid[0][j] = 1.0;
    }
    upc_forall( j=0; j<N+2; j++; j)
    {
        ptr[j] = &grid[j][0];
        new_ptr[j] = &new_grid[j][0];
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

    /* Set the precision wanted */
    epsilon  = 0.0001;
    finished = 0;
    nr_iter = 0;

    /* and start the timed section */
    gettimeofday( &ts_st, NULL );

    do
    {
        dTmax[MYTHREAD] = 0.0; 
        upc_forall( i=1; i<N+1; i++; i)
        {
            for( j=1; j<N+1; j++ )
            {
                T = 0.25 *
                    (ptr[i+1][j] + ptr[i-1][j] +
                     ptr[i][j-1] + ptr[i][j+1]); /* stencil */
                dT = T - ptr[i][j]; /* local variation */
                new_ptr[i][j] = T;
                if( dTmax[MYTHREAD] < fabs(dT) )
                    dTmax[MYTHREAD] = fabs(dT); /* max variation in this iteration */
            }
        }
        upc_barrier;
        upc_all_reduceD( &diffmax, dTmax, UPC_MAX,THREADS,1,NULL,UPC_OUT_ALLSYNC);
        if( diffmax < epsilon ) /* is the precision reached good enough ? */
            finished = 1;
        else
        {
            double shared *tmp;
            upc_forall( k=0; k<N+2; k++; k)      /* not yet ... Need to prepare */
            {                                    /* the next iteration */
                tmp = ptr[k];
                ptr[k] = new_ptr[k];
                new_ptr[k] = tmp;
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
