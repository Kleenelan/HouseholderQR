#include "utils.h"

void init_matrix(int M, int N, float* A, int lda, int seed)
{
    srand(seed);

    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){

            int r;
            r = rand();
            A[i + j*lda] = (((r>(RAND_MAX/3)) ? 3.0f : -3.0f)*rand())/RAND_MAX;
        }
    }
}

void print_matrix(int M, int N, float* A, int lda)
{
    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            printf("%7.4f, ", A[i + j*lda]);
        }
        printf("\n");
    }
}

void print_vector(int N, float* A)
{
    print_matrix(1, N, A, 1);
}