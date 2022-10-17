#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define N 4

double grid[N+2][N+2], new_grid[N+2][N+2];

void initialize(void)
{
    int j;

    /* Heat one side of the solid */
    for( j=1; j<N+1; j++ )
    {
        grid[0][j] = 1.0;
        new_grid[0][j] = 1.0;
    }
}

int main(void)
{
    struct timeval ts_st, ts_end;
    double dTmax, dT, epsilon, time;
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
        dTmax = 0.0; 
        for( i=1; i<N+1; i++ )
        {
            for( j=1; j<N+1; j++ )
            {
                T = 0.25 *
                    (grid[i+1][j] + grid[i-1][j] +
                     grid[i][j-1] + grid[i][j+1]); /* stencil */
                dT = T - grid[i][j]; /* local variation */
                new_grid[i][j] = T;
                if( dTmax < fabs(dT) )
                    dTmax = fabs(dT); /* max variation in this iteration */
            }
        }

        printf("nr_iter = %d, diffmax = %f", nr_iter, dTmax);
        print_grid();

        if( dTmax < epsilon ) /* is the precision reached good enough ? */
            finished = 1;
        else
        {
            for( k=0; k<N+2; k++ )      /* not yet ... Need to prepare */
                for( l=0; l<N+2; l++ )    /* ourselves for doing a new */
                    grid[k][l] = new_grid[k][l]; /* iteration */
        }
        nr_iter++;


    } while( finished == 0 );

    gettimeofday( &ts_end, NULL ); /* end the timed section */

    /* compute the execution time */
    time = ts_end.tv_sec + (ts_end.tv_usec / 1000000.0);
    time -= ts_st.tv_sec + (ts_st.tv_usec / 1000000.0);

    printf("%d iterations in %.3lf sec\n", nr_iter, time);

    print_grid();

    return 0;
}

int print_grid()
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

    return 0;

}
