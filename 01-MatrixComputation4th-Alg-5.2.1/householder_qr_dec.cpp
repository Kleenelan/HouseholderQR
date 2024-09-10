#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "householder_vector.h"

//
void vector_mul_matrix(int m, int n, float* vt, float* A, int lda, float* C)// C(1, n) = vt(1, m) * A(m, n)
{
    for(int j=0; j<n; j++){
        float sigma = 0.0f;
        for(int k=0; k<m; k++){
            sigma += vt[k] * A[k + j*lda];
        }
        C[j] = sigma;
    }
}

// A(m, n) = A - b * v(m) * C(n)
void rank_one_update(int m, int n, float beta, float* v, float* C, float* A, int lda)
{
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            A[i + j*lda] -= beta * v[i] * C[j];
        }
    }
}

void store_householder_vector(float* A, float* v, int len)
{
    for(int i=0; i<len; i++){
        A[i] = v[i];
    }
}

// (Householder QR)
void householder_qr(int m, int n, float* A, int lda, float* tau)
{
    float* v = nullptr;
    float beta = 0;
    float* C = nullptr;

    C = (float*)malloc(m*sizeof(float));
    v = (float*)malloc(m*sizeof(float));

    for(int j=0; j<n; j++){
        //step1, [v, beta] = house(A(j:m, j))
        //void house(int m, float* x, float* v, float& beta);
        house(m-j, A+j+j*lda, v, beta);
        tau[j] = beta;
        //step2, A(j:m, j:n) = (I - beta*v*v^T)A(j:m, j:n)//m 可能挺大，n<=32 or 64;
        // A = A - b*v*(v^t * A) = A - b*v*C;  (j:m, j:n),   v(j:m),   (v^t * A)(j:n)

          //step2.1 C(j:n) = (v^t)(1, j:m) * A(j:m, j:n)
        vector_mul_matrix(m-j, n-j, v, A+j+j*lda, lda, C);

          //step2.2 A(j:m, j:n) = A -  b*v*C;
        rank_one_update(m-j, n-j, beta, v, C, A+j+j*lda, lda);
        //step3, if(j<m) A(j+1 : m, j) = v(2 : m-j+1)
        store_householder_vector(A+(j+1 + j*lda), v+1, m-j-1);
    }
}



//#ifdef BUILD_MAIN
#if 1
int main()
{
    int m = 4;
    int n = 4;
    int lda = m;
    int ldb = lda;
    float* A = nullptr;
    float* B = nullptr;

    A = (float*)malloc(lda*n*sizeof(float));
    B = (float*)malloc(ldb*n*sizeof(float));

    init_matrix(m, n, A, lda, 2024); printf("A =[ ...\n");   print_matrix(m, n, A, lda);
    memcpy(B, A, lda*n*sizeof(float));

    float* tau = nullptr;
    tau = (float*)malloc(n*sizeof(float));

    householder_qr(m, n, A, lda, tau);
    printf("R+tau =\n");
    print_matrix(m, n, A, lda);

    printf("\ntau = \n");print_vector(n, tau);

    return 0;
}



#endif