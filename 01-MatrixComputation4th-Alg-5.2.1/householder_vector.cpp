#include <stdlib.h>


#include "utils.h"

//#define BUILD_MAIN
void v_is_1_x_2_m(int M, float* x, float* v)
{
    v[0] = 1.0f;
    for(int i=1; i<M; i++){
        v[i] = x[i];
    }
}

void dot_vector(int n, float* x, float* y, float& sigma)
{
    sigma = 0.0f;

    for(int i=0; i<n; i++)
    {
        sigma += x[i]*y[i];
    }
}

void vector_div_scalar(int M, float* v, float alpha)
{
    for(int i=0; i<M; i++){
        v[i] /= alpha;
    }
}

void house(int m, float* x, float* v, float& beta)
{
    float sigma;
    dot_vector(m-1, x+1, x+1, sigma);//    printf("in house() sigma = %7.3f\n", sigma);
    v_is_1_x_2_m(m, x, v);// v= ( 1 x(2 : m)^t )^t

    if(sigma==0.0f && x[0]>=0.0f)
    {
        beta = 0.0f;
    }
    else if(sigma==0.0f && x[0]<0.0f)
    {
        beta = 2.0f;
    }
    else
    {
        float miu;

        miu = sqrt(x[0]*x[0] + sigma);
        if(x[0]<= 0.0f){
            v[0] = x[0] - miu;
        }
        else{
            v[0] = -sigma/(x[0]+miu);
        }

        beta = 2.0f*v[0]*v[0]/(sigma + v[0]*v[0]);
        vector_div_scalar(m, v, v[0]);
    }
}

#ifdef BUILD_MAIN_VEC
//  #if 0
int main()
{
    int m = 7;
    float* x = nullptr;
//    printf("RAND_MAX = %d\n", RAND_MAX);

    x = (float*)malloc(m*sizeof(float));
    init_matrix(m, 1, x, m, 2024);
    print_vector(m, x);

    float sigma = 0.0;
    float* v = nullptr;

    v = (float*)malloc(m*sizeof(float));
    float beta;
    house(m, x, v, beta); printf("v =\n"); print_vector(m, v); printf("\nbeta = %7.4f\n", beta);

}

#endif
