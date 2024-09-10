#pragma once
#include <stdio.h>
#include <math.h>

void init_matrix(int M, int N, float* A, int lda, int seed);
void print_matrix(int M, int N, float* A, int lda);
void print_vector(int M, float* A);

