#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

long double w = 0.8;

void            relax_method( long double **A , long double *f , long double eps , int n );
void            one_iteraton( long double **A , long double *f , long double *x_k , long double *x_k_1 , int n );
long double     accuracy_of_solve_KOSHI( long double *x_k , long double *x_k_1 , int n );
long double     accuracy_of_solve_MATRIX( long double **A , long double *f , long double *x_k , int n );
int             read_matrix_system( const char *file_name , long double ***MATRIX , long double **STB );
void            write_matrix( const char *file_name , long double **A , long double *f , int n );
void            write_solution( const char *file_name , long double *f , int n );
void            make_matrix( const char *file_name , int mode );
void            free_mem( long double **A , long double *f , int n );

int
main( int argc , char **argv )
{
    //______________________
    //its only for use this packet of meth
    if ( argc == 2 ) {
        printf( "arg1 : input file\n" );
        printf( "arg2 : output file\n" );
        printf( "arg3 : accuracy\n" );
        printf( "arg4 : matrix in file <1> , matrix by fomula <2>\n" );
        return 0;
    }
    long double eps;
    sscanf( argv[3] , "%Lf" , &eps );
    int task;
    sscanf( argv[4] , "%d" , &task );
    unlink( argv[2] );
    if ( task == 1 ) {
        long double **A = NULL , *f = NULL;
        int size = read_matrix_system( argv[1] , &A , &f );
        relax_method( A , f , eps , size );
        write_solution( argv[2] , f , size );
        free_mem( A , f , size );
        return 0;
    }
    if ( task == 2 ) {
        long double **A = NULL , *f = NULL;
        make_matrix( argv[1] , 0 );
        int size = read_matrix_system( argv[1] , &A , &f );
        relax_method( A , f , eps , size );
        write_solution( argv[2] , f , size );
        free_mem( A , f , size );
        return 0;
    }
    if ( task == 3 ) {
        long double **A = NULL , *f = NULL;
        w = 0.05;
        make_matrix( argv[1] , 0 );
        for ( ; w < 1.97 ; w += 0.05 ) {
            int size = read_matrix_system( argv[1] , &A , &f );
            printf( "W = %Lf " , w );
            relax_method( A , f , eps , size );
            free_mem( A , f , size );
            A = NULL ; f = NULL;
        }
    }
    return 0;
}

void
relax_method( long double **A , long double *f , long double eps , int n )
{
    long double *x_k = (long double *) malloc ( n * sizeof( *x_k ) );
    long double *x_k_1 = (long double *) malloc ( n * sizeof( *x_k_1 ) );
    for ( int i = 0 ; i < n ; i++ ) { //initialized start value
        x_k[i] = 1;
    }
    //______________________
    //make next iteration of algo
    int k = 0;
    while ( 1 ) {
        if ( accuracy_of_solve_MATRIX( A , f , x_k , n ) < eps ) { // accuracy_of_solve_KOSHI( x_k , x_k_1 , n ) < eps
            break;
        }
        one_iteraton( A , f , x_k , x_k_1 , n );
        k++;
    }
    //======================
    //______________________
    //copy of solution in f
    for ( int i = 0 ; i < n ; i++ ) {
        f[i] = x_k[i];
    }
    //======================
    printf( "COUNT_IT :: %d\n" , k );
    free( x_k );
    free( x_k_1 );
}

void
one_iteraton( long double **A , long double *f , long double *x_k , long double *x_k_1 , int n )
{
    //______________________
    //iteration
    long double q = 1 - ( 1 / w ); //its koef of k_i
    for ( int i = 0 ; i < n ; i++ ) {
        x_k_1[i] = f[i];
        for ( int j = 0 ; j < i ; j++ ) {
            x_k_1[i] -= A[i][j] * x_k_1[j];
        }
        x_k_1[i] -= q * A[i][i] * x_k[i];
        for ( int j = i + 1 ; j < n ; j++ ) {
            x_k_1[i] -= A[i][j] * x_k[j];
        }
        x_k_1[i] *= w;
        x_k_1[i] /= A[i][i];
    }
    for ( int i = 0 ; i < n ; i++ ) {
        long double dop;
        dop = x_k[i];
        x_k[i] = x_k_1[i];
        x_k_1[i] = dop;
    }
    //really now i+i is x_k and i is x_k_1
    //======================
}

long double
accuracy_of_solve_KOSHI( long double *x_k , long double *x_k_1 , int n )
{
    //______________________
    //calc norma of solve
    long double norma = 0;
    for ( int i = 0 ; i < n ; i++ ) {
        norma += fabsl( x_k[i] - x_k_1[i] );
    }
    //======================
    return norma;
}

long double
accuracy_of_solve_MATRIX( long double **A , long double *f , long double *x_k , int n )
{
    //______________________
    //calc norma of solve
    long double *Ax = (long double *) calloc( n , sizeof( *Ax ) );
    for ( int i = 0 ; i < n ; i++ ) {
        for ( int j = 0 ; j < n ; j++ ) {
            Ax[i] += A[i][j] * x_k[j];
        }
    }
    long double norma = 0;
    for ( int i = 0 ; i < n ; i++ ) {
        norma += fabsl( f[i] - Ax[i] );
    }
    free( Ax );
    //======================
    return norma;
}

int
read_matrix_system( const char *file_name , long double ***MATRIX , long double **STB )
{
    //______________________
    //allocate size of matrix
    long double **A = *MATRIX;
    long double *f = *STB;
    FILE *in = fopen( file_name , "r" );
    int n;
    fscanf( in , "%d" , &n );
    // printf( "SIZE :: %d\n" , n );
    //======================
    //______________________
    //get memory
    A = (long double **) malloc( n * sizeof( *A ) );
    for ( int i = 0 ; i < n ; i++ ) {
        A[i] = (long double *) malloc ( n * sizeof( **A ) );
    }
    f = (long double *) malloc ( n * sizeof( *f ) );
    //======================
    //______________________
    //read matrix from file with file_name
    for ( int i = 0 ; i < n ; i++ ) {
        for ( int j = 0 ; j < n ; j++ ) {
            fscanf( in , "%Lf" , &A[i][j] );
        }
        fscanf( in , "%Lf" , &f[i] );
    }
    fclose( in );
    //======================
    // printf( "MATRIX OF SYSTEM WAS READED\n" );
    *MATRIX = A;
    *STB = f;
    return n;
}

void
write_matrix( const char *file_name , long double **A , long double *f , int n )
{
    //______________________
    //write matrix in file with file_name
    FILE *out = fopen( file_name , "a" );
    fprintf( out , "NEW MATRIX\n" );
    for ( int i = 0 ; i < n ; i++ ) {
        for ( int j = 0 ; j < n ; j++ ) {
            fprintf( out , "%Lf " , A[i][j] );
        }
        if ( f ) {
            fprintf( out , "| %Lf\n" , f[i] );
        } else {
            fprintf( out , "\n" );
        }
    }
    fprintf( out , "END MATRIX\n" );
    fclose( out );
    printf( "MATRIX WRITE IN FILE : %s\n" , file_name );
    //======================
}

void
write_solution( const char *file_name , long double *f , int n )
{
    //______________________
    //write solution only in file
    FILE *out = fopen( file_name , "a" );
    fprintf( out , "SOLUTION\n");
    for( int i = 0 ; i < n ; i++ ) {
        fprintf( out , "x_%d = %0.15Lf\n" , i + 1 , f[i] );
    }
    fclose( out );
    //======================
}

void
make_matrix( const char *file_name , int mode ) //gen matrix: mode 1 - only matrix ; 0 - full system
{
    //______________________
    //gen matrix
    int n = 30 , m = 20;
    FILE *out = fopen( file_name , "w" );
    fprintf( out , "%d\n" , n );
    for ( int i = 1 ; i <= n ; i++ ) {
        for ( int j = 1 ; j <= n ; j++ ) {
            if ( i != j ) {
                fprintf( out , "%Lf " , (long double) (i + j) / (long double) (n + m) );
            } else {
                fprintf( out , "%Lf " , (long double) n + ( m * m ) + ( (long double) j / m ) + ( (long double) i / n ) );
            }
        }
        if ( mode == 0 ) {
            fprintf( out , "%Lf\n" , (long double) m * i + n );
        } else {
            fprintf( out , "\n" );
        }
    }
    fclose( out );
    //======================
}

void
free_mem( long double **A , long double *f , int n )
{
    //______________________
    //free memory after use matrix
    for( int i = 0 ; i < n ; i++ ) {
        free( A[i] );
    }
    free( A );
    free( f );
    //======================
}
