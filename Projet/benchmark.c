/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#include "benchmark.h"

void benchmark_mul_naive_vs_strassen(u32 PRIME, u32 N, unsigned int ITE){
    printf("--- Multiplication : Naive vs Strassen ---\n");
    unsigned i, j;

    double time_naive = 0;
    double time_strassen = 0;
    double tic, toc, tic1, toc1;

    Matrix *A = init_matrix(PRIME, N);
    Matrix *B = init_matrix(PRIME, N);

    Matrix *res_stra;
    Matrix *res_naive;
    int k = 0;
    for(k = 0 ; k < ITE ; k++){
        for(i = 0 ; i < A->n ; i++){
            for(j = 0 ; j < A->n ; j++){
                A->m[i * A->n + j] = rand() % PRIME;
                B->m[i * B->n + j] = rand() % PRIME;
            }
        }
        // printf("--- STRASSEN ---\n");
        tic1 = wtime();
        res_stra = strassen_mul(A, B);
        // print_matrix(strassen_mul(A, B));
        toc1 = wtime();
        free_matrix(res_stra);
        time_strassen += toc1 - tic1;
        
        // printf("--- MUL ---\n");
        tic = wtime();
        res_naive = matrix_mul(A, B);
        // print_matrix(matrix_mul(A, B));
        toc = wtime();
        free_matrix(res_naive);
        time_naive += toc - tic;
    }
    free_matrix(A);
    free_matrix(B);

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- MATRIX MUL   = %f ---\n", time_naive);
    printf("--- STRASSEN MUL = %f ---\n", time_strassen);
    printf("##### Temps moyen #####\n");
    printf("--- MATRIX MUL   = %f ---\n", time_naive / ITE);
    printf("--- STRASSEN MUL = %f ---\n", time_strassen / ITE);

    printf("%s\n", (time_naive-time_strassen >0)? "Strassen meilleur" : "Naive meilleur");
}

void benchmark_inversion_pluq_vs_strassen(u32 PRIME, u32 N, unsigned int ITE){
    printf("--- Matrice inversion : PLUQ vs Strassen\n");
    unsigned i, j;

    double time_pluq = 0;
    double time_strassen = 0;
    double tic, toc, tic1, toc1;

    Matrix *A = init_matrix(PRIME, N);

    Matrix *res_stra;
    Matrix *res_pluq;

    int k = 0;
    for(k = 0 ; k < ITE ; k++){
        for(i = 0 ; i < A->n ; i++)
            for(j = 0 ; j < A->n ; j++)
                A->m[i * A->n + j] = (rand() % PRIME)+1;
            
        // printf("--- STRASSEN ---\n");
        tic1 = wtime();
        res_stra = strassen_inversion(A);
        // print_matrix(res_stra);
        toc1 = wtime();
        free_matrix(res_stra);
        time_strassen += toc1 - tic1;
        
        // printf("--- PLUQ ---\n");
        tic = wtime();
        res_pluq = pluq_inversion(A);
        // print_matrix(res_pluq);
        toc = wtime();
        free_matrix(res_pluq);
        time_pluq += toc - tic;
    }
    free_matrix(A);

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- INVERSION PLUQ     = %f ---\n", time_pluq);
    printf("--- INVERSION STRASSEN = %f ---\n", time_strassen);
    printf("##### Temps moyen #####\n");
    printf("--- INVERSION PLUQ     = %f ---\n", time_pluq / ITE);
    printf("--- INVERSION STRASSEN = %f ---\n", time_strassen / ITE);

    printf("%s\n", (time_pluq-time_strassen >0)? "Strassen meilleur" : "PLUQ meilleur");
}

/* Comparaison de l'inversion de Strassen avec la multiplication naive, strassen et optimise*/
void benchmark_strassen_inversion_StrassenMul_vs_naive_vs_optimized(u32 PRIME, u32 N, unsigned int ITE){
    printf("--- Matrice Strassen inversion : Strassen vs Strassen Multiplication vs Optimized (Mix) ---\n");
    unsigned i, j;

    double time_mul = 0;
    double time_naive = 0;
    double time_opti = 0;
    double tic, toc, tic1, toc1, tic2, toc2;

    Matrix *A = init_matrix(PRIME, N);

    Matrix *res_mul;
    Matrix *res_naive;
    Matrix *res_opti;

    int k = 0;
    for(k = 0 ; k < ITE ; k++){
        for(i = 0 ; i < A->n ; i++)
            for(j = 0 ; j < A->n ; j++)
                A->m[i * A->n + j] = (rand() % PRIME)+1;
            
        // printf("--- STRASSEN NAIVE ---\n");
        tic = wtime();
        res_naive = strassen_inversion(A);
        toc = wtime();
        // print_matrix(res_naive);
        free_matrix(res_naive);
        time_naive += toc - tic;
        
        // printf("--- STRASSEN MUL ---\n");
        tic1 = wtime();
        res_mul = strassen_mul_inversion(A);
        toc1 = wtime();
        // print_matrix(res_mul);
        free_matrix(res_mul);
        time_mul += toc1 - tic1;

        // printf("--- STRASSEN OPTI ---\n");
        tic2 = wtime();
        res_opti = strassen_optimized_inversion(A);
        toc2 = wtime();
        // print_matrix(res_opti);
        free_matrix(res_opti);
        time_opti += toc2 - tic2;
    }
    free_matrix(A);

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- INVERSION NAIVE = %f ---\n", time_naive);
    printf("--- INVERSION MUL   = %f ---\n", time_mul);
    printf("--- INVERSION OPTI  = %f ---\n", time_opti);
    printf("##### Temps moyen #####\n");
    printf("--- INVERSION NAIVE = %f ---\n", time_naive / ITE);
    printf("--- INVERSION MUL   = %f ---\n", time_mul / ITE);
    printf("--- INVERSION OPTI  = %f ---\n", time_opti / ITE);
}

void benchmark_pluq(u32 PRIME, u32 N, unsigned int ITE){
    printf("--- PLUQ ---\n");
    unsigned i, j;

    double time_pluq = 0;
    double tic, toc;

    Matrix *A = init_matrix(PRIME, N);

    PLUQ *res_pluq;

    int k = 0;
    for(k = 0 ; k < ITE ; k++){
        for(i = 0 ; i < A->n ; i++)
            for(j = 0 ; j < A->n ; j++)
                A->m[i * A->n + j] = (rand() % PRIME)+1;
            
        // printf("--- PLUQ ---\n");
        tic = wtime();
        res_pluq = pluq(A);
        toc = wtime();
        free_pluq(res_pluq);
        time_pluq += toc - tic;
    }

    free_matrix(A);

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- PLUQ = %f ---\n", time_pluq);
    printf("##### Temps moyen #####\n");
    printf("--- PLUQ = %f ---\n", time_pluq / ITE);
}


