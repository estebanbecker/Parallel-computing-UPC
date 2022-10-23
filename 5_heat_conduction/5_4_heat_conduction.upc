#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <upc_relaxed.h>
#include <upc_collective.h> 

#define grid(i,j) sh_grid[((i * N) + j)/(N*priv_size)].chunk[((i * N) + j)%(N*priv_size)]
#define new_grid(i,j) sh_new_grid[((i * N) + j)/(N*priv_size)].chunk[((i * N) + j)%(N*priv_size)]


shared double dTmax[THREADS];
shared double diffmax;
int N;
int priv_size;

typedef struct chunk_s chunk_t;
struct chunk_s {
    shared [] double *chunk;
};

shared chunk_t sh_grid[THREADS];
shared chunk_t sh_new_grid[THREADS];
chunk_t tmp;

void initialize(void)
{
    int j;

    /* Heat one side of the solid */
    for(j=1;j<N+1; j++)
    {
        grid(0,j) = 1.0;
        new_grid(0,j) = 1.0;
    }
}

int main(int argc, char *argv[])
{
    struct timeval ts_st, ts_end;
    double dT, epsilon, time;
    int finished, i, j, k, l;
    double T;
    int nr_iter;

    if(argc != 2)
    {
        if(MYTHREAD == 0)
        {
            printf("Usage: %s <N>", argv[0]);
        }        
        exit(1);
    }
    N = atoi(argv[1]);
    if(N<1)
    {
        if(MYTHREAD == 0)
        {
            printf("N must be greater than 1");
        }
        exit(1);
    }
    priv_size = (N+2)/THREADS;

    /*Alocate the memory*/
    sh_grid[MYTHREAD].chunk = (shared[] double *) upc_alloc((N+2)*priv_size*sizeof(double));
    sh_new_grid[MYTHREAD].chunk = (shared[] double *) upc_alloc((N+2)*priv_size*sizeof(double));

    if(MYTHREAD == 0)
    {
        initialize();
    }
    upc_barrier;
    /* Set the precision wanted */
    epsilon  = 0.0001;
    finished = 0;
    nr_iter = 0;

    /*Alloc the private pointers*/
    double **ptr_priv = (double **) malloc(priv_size*sizeof(double*));
    double **new_ptr_priv = (double **) malloc(priv_size*sizeof(double*));

    for(i = 0; i < priv_size; i++)
    {
        ptr_priv[i] = (double *) &grid(i + MYTHREAD*priv_size,0);
        new_ptr_priv[i] = (double *) &new_grid(i + MYTHREAD*priv_size,0);
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
                    (ptr_priv[i + 1][j] + grid((MYTHREAD*priv_size)-1,j) +
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
                    (grid((MYTHREAD + 1)* priv_size,j) + ptr_priv[i-1][j] +
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
            tmp = sh_grid[MYTHREAD];
            sh_grid[MYTHREAD] = sh_new_grid[MYTHREAD];
            sh_new_grid[MYTHREAD] = tmp;

            for( l=0; l<priv_size; l++)
            {
                double * tmp_priv = ptr_priv[l];
                ptr_priv[l] = new_ptr_priv[l];
                new_ptr_priv[l] = tmp_priv;
            }
        }
        upc_barrier;
        nr_iter++;
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
    upc_barrier;
    if(MYTHREAD==0)
    {
        printf("Grid:\n");
        int i, j;
        for( i=0; i<N+2; i++ )
        {
            for( j=0; j<N+2; j++ )
                printf("%f ", new_grid(i,j));
                scanf("%d", &i);
            printf("\n");
        }
        
    }
    
    upc_barrier;
    return 0;

}