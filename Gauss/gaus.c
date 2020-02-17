#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

long double eps = 0.000000001;

void            simple_gaus_front( long double **A , long double *f , int n );
void            gaus_reverse( long double **A , long double *f , int n );
void            mainelem_gaus( long double **A , long double *f , int n );
long double **  rev_matrix_gaus( long double **A , int n , long double *det );
long double     cond_num( long double **A , long double **B , int n );
void            swap_line( long double **A , long double *f , int f_1 , int f_2 , int n );
void            swap_stb( long double **A , int f_1 , int f_2 , int n );
int             read_matrix_system( const char *file_name , long double ***MATRIX , long double **STB );
int             read_matrix( const char *file_name , long double ***MATRIX );
void            write_matrix( const char *file_name , long double **A , long double *f , int n );
void            write_solution( const char *file_name , long double *f , int n );
long double **  copy_of_matrix( long double **res , long double **sen , int n );
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
        printf( "arg3 : task [1 , 2 , 3 , 4 , 5 , 6]\n" );
        printf( "1 - solve by simple_gaus\n" );
        printf( "2 - solve by mainelem_gaus\n" );
        printf( "3 - find A^-1\n" );
        printf( "4 - find ObNum\n" );
        printf( "5 - find det\n" );
        printf( "6 - make matrix by formula and write it to file arg1 , and arg4 is task [1 - 5]" );
        return 0;
    }
    int task = 0;
    sscanf( argv[3] , "%d" , &task );
    if ( task == 6 ) {
        sscanf( argv[4] , "%d" , &task );
        if ( task >= 3 && task <= 5 ) {
            make_matrix( argv[1] , 1 );
        } else {
            make_matrix( argv[1] , 0 );
        }
    }
    unlink( argv[2] );
    if ( task == 1 ) {
        long double **A = NULL , *f = NULL;
        int size = read_matrix_system( argv[1] , &A , &f );
        simple_gaus_front( A , f , size );
        gaus_reverse( A , f , size );
        write_solution( argv[2] , f , size );
        free_mem( A , f , size );
        return 0;
    }
    if ( task == 2 ) {
        long double **A = NULL , *f = NULL;
        int size = read_matrix_system( argv[1] , &A , &f );
        mainelem_gaus( A , f , size );
        write_solution( argv[2] , f , size );
        free_mem( A , f , size );
        return 0;
    }
    if ( task == 3 ) {
        long double **A = NULL;
        long double det;
        int size = read_matrix( argv[1] , &A );
        long double **B = rev_matrix_gaus( A , size , &det );
        write_matrix( argv[2] , B , NULL , size );
        free_mem( A , NULL , size );
        free_mem( B , NULL , size );
        return 0;
    }
    if ( task == 4 ) {
        long double **A = NULL;
        long double det;
        int size = read_matrix( argv[1] , &A );
        long double **B = rev_matrix_gaus( A , size , &det );
        long double ob = cond_num( A , B , size );
        FILE *out = fopen( argv[2] , "a" );
        fprintf( out , "ObNum: %Lf\n" , ob );
        fclose( out );
        free_mem( A , NULL , size );
        free_mem( B , NULL , size );
        return 0;
    }
    if ( task == 5 ) {
        long double **A = NULL;
        long double det;
        int size = read_matrix( argv[1] , &A );
        long double **B = rev_matrix_gaus( A , size , &det );
        FILE *out = fopen( argv[2] , "a" );
        fprintf( out , "Det: %Lf\n" , det );
        fclose( out );
        free_mem( A , NULL , size );
        free_mem( B , NULL , size );
        return 0;
    }
    //======================
    return 0;
}

void
simple_gaus_front( long double **A , long double *f , int n )
{
    printf( "GAUS_FRONT WAS STARTED\n" );
    //______________________
    //the list to the top go
    for ( int i = 0 ; i < n ; i++ ) { //string iterator
        //______________________
        //find non zero string and mode matrix
        int num_non_zero = i;
        for ( ; num_non_zero < n ; num_non_zero++ ) { //find not zero elem in stb
            if ( fabsl( A[num_non_zero][i] ) > eps ) {
                break;
            }
        }
        if ( i != num_non_zero ) { //if not zero elem in set position
            swap_line( A , f , i , num_non_zero , n );
        }
        //======================
        //______________________
        //devide string
        for ( int j = i + 1 ; j < n ; j++ ) { //stb iterator
            A[i][j] /= A[i][i]; //take main num
        }
        f[i] /= A[i][i];
        A[i][i] = 1;
        //======================
        //______________________
        //reduse matrix with round result
        for ( int k = i + 1 ; k < n ; k++ ) { //dec of the string to zero start number
            for ( int p = i + 1 ; p < n ; p++ ) {
                A[k][p] -= A[i][p] * A[k][i];
            }
            f[k] -= f[i] * A[k][i];
            A[k][i] = 0;
        }
        //======================
    }
    //======================
    printf( "GAUS_FRONT WAS FINISHED\n" );
}

void
gaus_reverse( long double **A , long double *f , int n )
{
    printf( "GAUS_REVERSE WAS STARTED\n" );
    //______________________
    //the list to go the end
    for ( int i = n - 1 ; i >= 0 ; i-- ) { //set in the last point in matrix
        for ( int j = i - 1 ; j >= 0 ; j-- ) {
            f[j] -= A[j][i] * f[i];
            A[j][i] = 0;
        }
    }
    //======================
    printf( "GAUS_REVERSE WAS FINISHED\n" );
}

void
mainelem_gaus( long double **A , long double *f , int n )
{
    printf( "MAINELEM_GAUS WAS STARTED\n" );
    //______________________
    //map of stb location for allocate start position
    int *map_of_stb_location = (int *) malloc( n * sizeof( *map_of_stb_location ) );
    for ( int i = 0 ; i < n ; i++ ) {
        map_of_stb_location[i] = i;
    }
    //======================
    //______________________
    //only gaus
    for ( int i = 0 ; i < n ; i++ ) { //string iterator
        //______________________
        //find max element
        int num_of_max = i; //num of max stb
        long double max = fabsl( A[i][i] );
        for ( int j = i + 1 ; j < n ; j++ ) { //find max elem in string
            if ( fabsl( A[i][j] ) > max ) {
                max = fabsl( A[i][j] );
                num_of_max = j;
            }
        }
        if ( i != num_of_max ) { //if not zero elem in set position
            swap_stb( A , i , num_of_max , n );
            int dop = map_of_stb_location[i];
            map_of_stb_location[i] = map_of_stb_location[num_of_max];
            map_of_stb_location[num_of_max] = dop;
        }
        //======================
        //______________________
        //devide string
        for ( int j = i + 1 ; j < n ; j++ ) { //stb iterator
            A[i][j] /= A[i][i]; //take main num
        }
        f[i] /= A[i][i];
        A[i][i] = 1;
        //======================
        //______________________
        //reduse matrix with round result
        for ( int k = i + 1 ; k < n ; k++ ) { //dec of the string to zero start number
            for ( int p = i + 1 ; p < n ; p++ ) {
                A[k][p] -= A[i][p] * A[k][i];
            }
            f[k] -= f[i] * A[k][i];
            A[k][i] = 0;
        }
        //======================
    }
    //======================
    //______________________
    //make finish deal
    gaus_reverse( A , f , n ); //start reverse gaus to make a solution
    long double *ask = (long double *) malloc( n * sizeof( *ask ) );
    for ( int i = 0 ; i < n ; i++ ) { //get start allocate base on start map
        ask[map_of_stb_location[i]] = f[i];
    }
    for ( int i = 0 ; i < n ; i++ ) { //copy the solution
        f[i] = ask[i];
    }
    free( ask ); //free mem
    free( map_of_stb_location ); //free mem
    //======================
    printf( "MAINELEM_GAUS WAS FINISHED\n" );
}

long double **
rev_matrix_gaus( long double **A , int n , long double *det )
{
    printf( "REV_MATRIX WAS STARTED\n" );
    *det = 1; //its only fo calc determ
    long double **copy_of_A = copy_of_matrix( A , NULL , n );
    //______________________
    //map of stb location for allocate start position
    int *map_of_stb_location = (int *) malloc( n * sizeof( *map_of_stb_location ) );
    for ( int i = 0 ; i < n ; i++ ) {
        map_of_stb_location[i] = i;
    }
    long double **B = (long double **) calloc( n , sizeof( *B ) );
    for( int i = 0 ; i < n ; i++ ) {
        B[i] = (long double *) calloc( n , sizeof( **B ) );
        B[i][i] = 1;
    }
    //======================
    //______________________
    //only gaus
    for ( int i = 0 ; i < n ; i++ ) { //string iterator
        //______________________
        //find max element
        int num_of_max = i; //num of max stb
        long double max = fabsl( A[i][i] );
        for ( int j = i + 1 ; j < n ; j++ ) { //find max elem in string
            if ( fabsl( A[i][j] ) > max ) {
                max = fabsl( A[i][j] );
                num_of_max = j;
            }
        }
        if ( i != num_of_max ) { //if not zero elem in set position
            *det *= -1; //its only fo calc determ
            swap_stb( A , i , num_of_max , n );
            int dop = map_of_stb_location[i];
            map_of_stb_location[i] = map_of_stb_location[num_of_max];
            map_of_stb_location[num_of_max] = dop;
        }
        //======================
        //______________________
        //devide string
        *det *= A[i][i]; //its only fo calc determ
        for ( int j = i + 1 ; j < n ; j++ ) { //stb iterator
            A[i][j] /= A[i][i]; //take main num
        }
        for ( int j = 0 ; j < n ; j++ ) {
            B[i][j] /= A[i][i];
        }
        A[i][i] = 1;
        //======================
        //______________________
        //reduse matrix with round result
        for ( int k = i + 1 ; k < n ; k++ ) { //dec of the string to zero start number
            for ( int p = i + 1 ; p < n ; p++ ) {
                A[k][p] -= A[i][p] * A[k][i];
            }
            for ( int p = 0 ; p < n ; p++ ) {
                B[k][p] -= B[i][p] * A[k][i];
            }
            A[k][i] = 0;
        }
        //======================
    }
    //======================
    //______________________
    //the list to go the end reverse go
    for ( int i = n - 1 ; i >= 0 ; i-- ) { //set in the last point in matrix
        for ( int j = i - 1 ; j >= 0 ; j-- ) {
            for ( int p = 0 ; p < n ; p++ ) {
                B[j][p] -= A[j][i] * B[i][p];
            }
            A[j][i] = 0;
        }
    }
    //======================
    //______________________
    //make finish deal
    long double **ask = (long double **) malloc( n * sizeof( *ask ) ); //swap string to male a result
    for ( int i = 0 ; i < n ; i++ ) { //get start allocate base on start map
        ask[map_of_stb_location[i]] = B[i];
    }
    for ( int i = 0 ; i < n ; i++ ) { //copy the solution
        B[i] = ask[i];
    }
    free( ask ); //free mem
    free( map_of_stb_location ); //free mem
    //======================
    A = copy_of_matrix( copy_of_A , A , n );
    free_mem( copy_of_A , NULL , n );
    printf( "REV_MATRIX WAS FINISHED\n" );
    return B;
}

long double
cond_num( long double **A , long double **B , int n )
{
    long double ob = 0;
    long double *mas_A = (long double *) calloc( n , sizeof( *mas_A ) );//mass of string sum in matrix A
    long double *mas_B = (long double *) calloc( n , sizeof( *mas_B ) );//mass of string sum in matrix B = A^-1
    //______________________
    //only find string sum of both matrix
    for( int i = 0 ; i < n ; i++ ) {
        for( int j = 0 ; j < n ; j++ ) {
            mas_A[i] += fabsl( A[i][j] );
            mas_B[i] += fabsl( B[i][j] );
        }
    }
    //======================
    //______________________
    //only find max sum of string of both matrix
    long double norma_A = mas_A[0]; //inf norma
    long double norma_B = mas_B[0]; //inf norma
    for ( int i = 1 ; i < n ; i++ ) {
        if ( mas_A[i] > norma_A ) {
            norma_A = mas_A[i];
        }
        if ( mas_B[i] > norma_B ) {
            norma_B = mas_B[i];
        }
    }
    //======================
    ob = norma_A * norma_B;
    free( mas_A );
    free( mas_B );
    return ob;
}

void
swap_line( long double **A , long double *f , int f_1 , int f_2 , int n ) //swap system line
{
    long double dop;
    for ( int p = f_1 ; p < n ; p++ ) { //swap string of matrix
        dop = A[f_1][p];
        A[f_1][p] = A[f_2][p];
        A[f_2][p] = dop;
    }
    if ( f ) {
        dop = f[f_1]; //swap string of f stb
        f[f_1] = f[f_2];
        f[f_2] = dop;
    }
}

void
swap_stb( long double **A , int f_1 , int f_2 , int n ) //swap A stb
{
    long double dop;
    for ( int p = 0 ; p < n ; p++ ) { //swap stb of matrix
        dop = A[p][f_1];
        A[p][f_1] = A[p][f_2];
        A[p][f_2] = dop;
    }
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
    printf( "SIZE :: %d\n" , n );
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
    printf( "MATRIX OF SYSTEM WAS READED\n" );
    *MATRIX = A;
    *STB = f;
    return n;
}

int
read_matrix( const char *file_name , long double ***MATRIX )
{
    //______________________
    //allocate size of matrix
    long double **A = *MATRIX;
    FILE *in = fopen( file_name , "r" );
    int n;
    fscanf( in , "%d" , &n );
    printf( "SIZE :: %d\n" , n );
    //======================
    //______________________
    //get memory
    A = (long double **) malloc( n * sizeof( *A ) );
    for ( int i = 0 ; i < n ; i++ ) {
        A[i] = (long double *) malloc ( n * sizeof( **A ) );
    }
    //======================
    //______________________
    //read matrix from file with file_name
    for ( int i = 0 ; i < n ; i++ ) {
        for ( int j = 0 ; j < n ; j++ ) {
            fscanf( in , "%Lf" , &A[i][j] );
        }
    }
    fclose( in );
    //======================
    printf( "MATRIX WAS READED\n" );
    *MATRIX = A;
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

long double **
copy_of_matrix( long double **A , long double **B , int n ) //its only copy matrix for save from algo
{
    //______________________
    //copy
    if ( !B ) {
        B = (long double **) calloc( n , sizeof( *B ) );
        for( int i = 0 ; i < n ; i++ ) {
            B[i] = (long double *) calloc( n , sizeof( **B ) );
        }
    }
    for ( int i = 0 ; i < n ; i++ ) {
        for ( int j = 0 ; j < n ; j++ ) {
            B[i][j] = A[i][j];
        }
    }
    //======================
    return B;
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
