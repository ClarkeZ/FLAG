/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#include "tests.h"

/* ---------- TEST BASE.C ---------- */

void test_add(u32 PRIME){
    printf("--- TEST ADD ---\n");
    u32 a = rand() % PRIME;
    u32 b = rand() % PRIME;
    u32 n = rand() % PRIME;

    printf("%u + %u = %u %% %u\n",a, b, add(a, b, n), n);
}

void test_sub(u32 PRIME){
    printf("--- TEST SUB ---\n");
    u32 a = rand() % PRIME;
    u32 b = rand() % PRIME;
    u32 n = rand() % PRIME;

    printf("%u - %u = %u %% %u\n",a, b, sub(a, b, n), n);
}

void test_mul(u32 PRIME){
    printf("--- TEST MUL ---\n");
    u32 a = rand() % PRIME;
    u32 b = rand() % PRIME;
    u32 n = rand() % PRIME;

    printf("%u * %u = %u %% %u\n",a, b, mul(a, b, n), n);
}

void test_mulmod(u32 PRIME){
    printf("--- TEST MULMOD ---\n");
    u32 a = rand() % PRIME;
    u32 b = rand() % PRIME;
    u32 n = rand() % PRIME;

    printf("%u * %u = %u %% %u\n",a, b, mulmod(a, b, n), n);
}

void test_invmod(u32 PRIME){
    printf("--- TEST INVMOD ---\n");
    u32 a = rand() % PRIME;
    u32 b = rand() % PRIME;

    printf("%u ^-1 [%u] = %u\n",a, b, invmod(a, b));
}

void test_exp_mod(u32 PRIME){
    printf("--- TEST EXP_MOD ---\n");
    u32 a = rand() % PRIME;
    u32 b = rand() % PRIME;
    u32 n = rand() % PRIME;

    printf("%u^%u mod %u = %u\n",a, b, n, exp_mod(a, b, n));
}



/* ---------- TEST MATRIX.C ---------- */

void test_init_matrix(u32 PRIME, u32 N){
    printf("--- TEST INIT_MATRIX ---\n");
    Matrix *M = init_matrix(PRIME, N);
    print_matrix(M);

    free_matrix(M);
}

void test_init_eye(u32 PRIME, u32 N){
    printf("--- TEST INIT_EYE ---\n");
    Matrix *M = init_eye(PRIME, N);
    print_matrix(M);

    free_matrix(M);
}

void test_copy_matrix(u32 PRIME, u32 N){
    printf("--- TEST COPY_MATRIX ---\n");
    Matrix *M = init_matrix(PRIME, N);
    Matrix *copy;
    double tic, toc, time;
    unsigned int i;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- M AVANT ---\n");
    print_matrix(M);

    tic = wtime();
    copy = copy_matrix(M);
    toc = wtime();

    printf("--- M APRES ---\n");
    print_matrix(M);
    printf("--- COPY APRES ---\n");
    print_matrix(copy);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
    free_matrix(copy);
}

void test_matrix_add(u32 PRIME, u32 N){
    printf("--- TEST MATRIX_ADD ---\n");
    Matrix *A = init_matrix(PRIME, N);
    Matrix *B = init_matrix(PRIME, N);
    Matrix *res;
    double tic, toc, time;;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = rand() % PRIME;
        B->m[i] = rand() % PRIME;
    }
    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- B ---\n");
    // print_matrix(B);

    tic = wtime();
    res = matrix_add(A, B);
    toc = wtime();

    // printf("--- A + B ---\n");
    // print_matrix(res);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(A);
    free_matrix(B);
    free_matrix(res);
}

void test_matrix_sub(u32 PRIME, u32 N){
    printf("--- TEST MATRIX_SUB ---\n");
    Matrix *A = init_matrix(PRIME, N);
    Matrix *B = init_matrix(PRIME, N);
    Matrix *res;
    double tic, toc, time;;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = rand() % PRIME;
        B->m[i] = rand() % PRIME;
    }
    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- B ---\n");
    // print_matrix(B);

    tic = wtime();
    res = matrix_sub(A, B);
    toc = wtime();

    // printf("--- A - B ---\n");
    // print_matrix(res);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(A);
    free_matrix(B);
    free_matrix(res);
}

void test_matrix_mul(u32 PRIME, u32 N){
    printf("--- TEST MATRIX_MUL ---\n");
    Matrix *A = init_matrix(PRIME, N);
    Matrix *B = init_matrix(PRIME, N);
    Matrix *res;
    double tic, toc, time;;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = rand() % PRIME;
        B->m[i] = rand() % PRIME;
    }
    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- B ---\n");
    // print_matrix(B);

    tic = wtime();
    res = matrix_mul(A, B);
    toc = wtime();

    // printf("--- A * B ---\n");
    // print_matrix(res);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(A);
    free_matrix(B);
    free_matrix(res);
}

void test_strassen_mul(u32 PRIME, u32 N){
    printf("--- TEST STRASSEN_MUL ---\n");
    Matrix *A = init_matrix(PRIME, N);
    Matrix *B = init_matrix(PRIME, N);
    Matrix *res;
    double tic, toc, time;;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = rand() % PRIME;
        B->m[i] = rand() % PRIME;
    }
    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- B ---\n");
    // print_matrix(B);

    tic = wtime();
    res = strassen_mul(A, B);
    toc = wtime();

    // printf("--- A * B ---\n");
    // print_matrix(res);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(A);
    free_matrix(B);
    free_matrix(res);
}

void test_matrix_mul_coef(u32 PRIME, u32 N){
    printf("--- TEST MATRIX_MUL_COEF ---\n");
    Matrix *A = init_matrix(PRIME, N);
    u32 c = rand() % PRIME;
    Matrix *res;
    double tic, toc, time;;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i)
        A->m[i] = rand() % PRIME;

    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- c ---\n");
    // printf("%u\n", c);

    tic = wtime();
    res = matrix_mul_coef(A, c);
    toc = wtime();

    // printf("--- A * c ---\n");
    // print_matrix(res);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(A);
    free_matrix(res);
}

void test_matrix_mul_optimized(u32 PRIME, u32 N){
    printf("--- TEST MATRIX_MUL_OPTIMIZED ---\n");
    Matrix *A = init_matrix(PRIME, N);
    Matrix *B = init_matrix(PRIME, N);
    Matrix *res;
    double tic, toc, time;;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = rand() % PRIME;
        B->m[i] = rand() % PRIME;
    }
    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- B ---\n");
    // print_matrix(B);

    tic = wtime();
    res = matrix_mul_optimized(A, B);
    toc = wtime();

    // printf("--- A * B ---\n");
    // print_matrix(res);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(A);
    free_matrix(B);
    free_matrix(res);
}

void test_transpose(u32 PRIME, u32 N){
    printf("--- TEST TRANSPOSE ---\n");
    Matrix *A = init_matrix(PRIME, N);
    Matrix *transp;
    double tic, toc, time;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i)
        A->m[i] = rand() % PRIME;

    printf("--- A ---\n");
    print_matrix(A);

    tic = wtime();
    transp = transpose(A);
    toc = wtime();

    printf("--- A^T ---\n");
    print_matrix(transp);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(transp);
    free_matrix(A);
}

/* ---------- TEST ALGO.C ---------- */

void test_swap_lines(u32 PRIME, u32 N){
    printf("--- TEST SWAP_LINES ---\n");
    Matrix *M = init_matrix(PRIME, N);
    int i;
    int l1 = 1;
    int l2 = 3;
    double tic, toc, time;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- AVANT ---\n");
    print_matrix(M);

    tic = wtime();
    swap_lines(M, l1, l2);
    toc = wtime();

    printf("--- APRES ---\n");
    print_matrix(M);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
}

void test_swap_columns(u32 PRIME, u32 N){
    printf("--- TEST SWAP_COLUMNS ---\n");
    Matrix *M = init_matrix(PRIME, N);
    int i;
    int c1 = 1;
    int c2 = 3;
    double tic, toc, time;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- AVANT ---\n");
    print_matrix(M);
    
    tic = wtime();
    swap_columns(M, c1, c2);
    toc = wtime();

    printf("--- APRES ---\n");
    print_matrix(M);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
}

void test_swap_pluq_pivot(u32 PRIME, u32 N){
    printf("--- TEST SWAP_PLUQ_PIVOT ---\n");
    PLUQ *pluqed = (PLUQ*) malloc(sizeof(PLUQ));
    pluqed->P = init_eye(PRIME, N);
    pluqed->L = init_eye(PRIME, N);
    pluqed->U = init_matrix(PRIME, N);
    pluqed->Q = init_eye(PRIME, N);

    int i;
    double tic, toc, time;

    for (i = 0 ; i < pluqed->U->n * pluqed->U->n ; ++i) 
        pluqed->U->m[i] = rand() % PRIME;

    printf("--- AVANT ---\n");
    printf("--- P ---\n");
    print_matrix(pluqed->P);
    printf("--- L ---\n");
    print_matrix(pluqed->L);
    printf("--- U ---\n");
    print_matrix(pluqed->U);
    printf("--- Q ---\n");
    print_matrix(pluqed->Q);

    tic = wtime();
    swap_pluq_pivot(pluqed, 0);
    toc = wtime();

    printf("--- APRES ---\n");
    printf("--- P ---\n");
    print_matrix(pluqed->P);
    printf("--- L ---\n");
    print_matrix(pluqed->L);
    printf("--- U ---\n");
    print_matrix(pluqed->U);
    printf("--- Q ---\n");
    print_matrix(pluqed->Q);

    time = toc - tic;
    printf("time = %f\n", time);

    free_pluq(pluqed);
}

void test_mul_line(u32 PRIME, u32 N){
    printf("--- TEST MUL_LINE ---\n");
    Matrix *M = init_matrix(PRIME, N);
    int i;
    int ligne = 2;
    int c = rand() % PRIME;
    double tic, toc, time;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- AVANT ---\n");
    print_matrix(M);
    printf("coefficient = %u\n", c);
    
    tic = wtime();
    mul_line(M, c, ligne);
    toc = wtime();

    printf("--- APRES ---\n");
    print_matrix(M);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
}

void test_multiply_line(u32 PRIME, u32 N){
    printf("--- TEST MULTIPLY_LINE ---\n");
    u32 *u = (u32 *) malloc(sizeof(u32) * N);
    u32 *res = (u32 *) malloc(sizeof(u32) * N);
    int i;
    int c = rand() % PRIME;
    double tic, toc, time;

    for (i = 0 ; i < N ; ++i) 
        u[i] = rand() % PRIME;

    printf("--- AVANT ---\nu = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",u[i]);
    printf("\ncoefficient = %u\n", c);
    
    tic = wtime();
    res = multiply_line(u, c, N, PRIME);
    toc = wtime();

    printf("--- APRES ---\nu = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",res[i]);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free(u);
    free(res);
}

void test_substract_lines(u32 PRIME, u32 N){
    printf("--- TEST SUBSTRACT_LINES ---\n");
    u32 *u = (u32 *) malloc(sizeof(u32) * N);
    u32 *v = (u32 *) malloc(sizeof(u32) * N);
    u32 *res = (u32 *) malloc(sizeof(u32) * N);
    int i;
    double tic, toc, time;

    for (i = 0 ; i < N ; ++i){
        u[i] = rand() % PRIME;
        v[i] = rand() % PRIME;
    }
    printf("--- AVANT ---\nu = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",u[i]);
    printf("\nv = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",v[i]);
    
    tic = wtime();
    res = substract_lines(u, v, N, PRIME);
    toc = wtime();

    printf("\n--- APRES ---\nu = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",res[i]);
    printf("\nv = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",v[i]);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free(u);
    free(v);
    free(res);
}

void test_sub_lines(u32 PRIME, u32 N){
    printf("--- TEST SUB_LINES ---\n");
    Matrix *M = init_matrix(PRIME, N);
    int i;
    int ligne1 = 2;
    int ligne2 = 3;
    double tic, toc, time;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- AVANT ---\n");
    print_matrix(M);
    
    tic = wtime();
    sub_lines(M, ligne1, ligne2);
    toc = wtime();

    printf("--- APRES ---\n");
    print_matrix(M);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
}

void test_zero_in_diagonal(u32 PRIME, u32 N){
    printf("--- TEST ZERO_IN_DIAGONAL ---\n");
    Matrix *M = init_matrix(PRIME, N);
    int i;
    double tic, toc, time;
    bool iszero;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- M ---\n");
    print_matrix(M);
    
    tic = wtime();
    iszero = zero_in_diagonal(M);
    toc = wtime();

    printf("%s\n", (iszero==true)?"True":"False");

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
}

void test_solv_upper_triangular_matrix(u32 PRIME, u32 N){
    printf("--- TEST SOLV_UPPER_TRIANGULAR_MATRIX ---\n");
    Matrix *M = init_matrix(PRIME, N);
    u32 *u = (u32 *) malloc(sizeof(u32) * N);
    u32 *res = (u32 *) malloc(sizeof(u32) * N);
    int i, j;
    double tic, toc, time;

    for (i = 0 ; i < N ; ++i) 
        for(j = i ; j < N ; ++j)
            M->m[i * N + j] = rand() % PRIME;

    for (i = 0 ; i < N ; ++i)
        u[i] = rand() % PRIME;
    printf("--- A ---\n");
    print_matrix(M);
    printf("b = ");
    for (i = 0 ; i < N ; ++i)
        printf("%u ",u[i]);
    
    tic = wtime();
    res = solv_upper_triangular_matrix(M, u);
    toc = wtime();

    printf("\nx = ");
    for (i = 0 ; i < N ; ++i)
        printf("%u ",res[i]);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free_matrix(M);
    free(u);
    free(res);
}

void test_solv_lower_triangular_matrix(u32 PRIME, u32 N){
    printf("--- TEST SOLV_LOWER_TRIANGULAR_MATRIX ---\n");
    Matrix *M = init_eye(PRIME, N);
    u32 *u = (u32 *) malloc(sizeof(u32) * N);
    u32 *res = (u32 *) malloc(sizeof(u32) * N);
    int i, j;
    double tic, toc, time;

    for (i = 0 ; i < N ; ++i) 
        for(j = 0 ; j < i ; ++j)
            M->m[i * N + j] = rand() % PRIME;

    for (i = 0 ; i < N ; ++i)
        u[i] = rand() % PRIME;
    printf("--- A ---\n");
    print_matrix(M);
    printf("b = ");
    for (i = 0 ; i < N ; ++i)
        printf("%u ",u[i]);
    
    tic = wtime();
    res = solv_lower_triangular_matrix(M, u);
    toc = wtime();

    printf("\nx = ");
    for (i = 0 ; i < N ; ++i)
        printf("%u ",res[i]);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free_matrix(M);
    free(u);
    free(res);
}

void test_inverse_upper_triangular_matrix(u32 PRIME, u32 N){
    printf("--- TEST INVERSE_UPPER_TRIANGULAR_MATRIX ---\n");
    Matrix *M = init_matrix(PRIME, N);
    Matrix *res;
    int i, j;
    double tic, toc, time;

    for (i = 0 ; i < N ; ++i) 
        for(j = i ; j < N ; ++j)
            M->m[i * N + j] = rand() % PRIME;

    printf("--- M ---\n");
    print_matrix(M);
    
    tic = wtime();
    res = inverse_upper_triangular_matrix(M);
    toc = wtime();

    printf("--- M^-1 ---\n");
    print_matrix(res);

    printf("--- VERIF ---\n");
    Matrix *tmp_strassen = strassen_mul(M, res);
    print_matrix(tmp_strassen);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free_matrix(tmp_strassen);
    free_matrix(M);
    free_matrix(res);
}

void test_inverse_lower_triangular_matrix(u32 PRIME, u32 N){
    printf("--- TEST INVERSE_LOWER_TRIANGULAR_MATRIX ---\n");
    Matrix *M = init_eye(PRIME, N);
    Matrix *res;
    int i, j;
    double tic, toc, time;

    for (i = 0 ; i < N ; ++i) 
        for(j = 0 ; j < i ; ++j)
            M->m[i * N + j] = rand() % PRIME;

    printf("--- M ---\n");
    print_matrix(M);
    
    tic = wtime();
    res = inverse_lower_triangular_matrix(M);
    toc = wtime();

    printf("--- M^-1 ---\n");
    print_matrix(res);

    printf("--- VERIF ---\n");
    Matrix *tmp_strassen = strassen_mul(M, res);
    print_matrix(tmp_strassen);
    free_matrix(tmp_strassen);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free_matrix(M);
    free_matrix(res);
}

void test_matrix_vector_product(u32 PRIME, u32 N){
    printf("--- TEST MATRIX_VECTOR_PRODUCT ---\n");
    Matrix *M = init_matrix(PRIME, N);
    u32 *u = (u32 *) malloc(sizeof(u32) * N);
    u32 *res = (u32 *) malloc(sizeof(u32) * N);
    int i;
    double tic, toc, time;

    for (i = 0 ; i < N * N ; ++i) 
        M->m[i] = rand() % PRIME;
    for (i = 0 ; i < N ; ++i)
        u[i] = rand() % PRIME;

    printf("--- AVANT ---\n");
    print_matrix(M);
    printf("u = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",u[i]);
    
    tic = wtime();
    res = matrix_vector_product(M, u);
    toc = wtime();

    printf("\n--- APRES ---\nres = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",res[i]);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free_matrix(M);
    free(u);
    free(res);
}

void test_pluq(u32 PRIME, u32 N){
    printf("--- TEST PLUQ ---\n");
    Matrix *M = init_matrix(PRIME, N);
    PLUQ *pluqed = (PLUQ *) malloc(sizeof(PLUQ));
    int i;
    double tic, toc, time;

    for (i = 0 ; i < N * N ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- AVANT ---\n");
    print_matrix(M);
    
    tic = wtime();
    pluqed = pluq(M);
    toc = wtime();

    printf("--- APRES ---\n");
    printf("--- P ---\n");
    print_matrix(pluqed->P);
    printf("--- L ---\n");
    print_matrix(pluqed->L);
    printf("--- U ---\n");
    print_matrix(pluqed->U);
    printf("--- Q ---\n");
    print_matrix(pluqed->Q);

    printf("--- VERIF = AVANT ---\n");
    Matrix *tmp_pluq1 = matrix_mul(pluqed->U,pluqed->Q);
    Matrix *tmp_pluq2 = matrix_mul(pluqed->L,tmp_pluq1);
    Matrix *tmp_pluq3 = matrix_mul(pluqed->P,tmp_pluq2);
    print_matrix(tmp_pluq3);

    free_matrix(tmp_pluq1);
    free_matrix(tmp_pluq2);
    free_matrix(tmp_pluq3);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free_matrix(M);
    free_pluq(pluqed);
}

void test_linear_solv(u32 PRIME, u32 N){
    printf("--- TEST LINEAR_SOLV ---\n");
    Matrix *M = init_matrix(PRIME, N);
    u32 *u = (u32 *) malloc(sizeof(u32) * N);
    u32 *res;
    int i;
    double tic, toc, time;

    for (i = 0 ; i < N * N ; ++i) 
        M->m[i] = rand() % PRIME;
    for (i = 0 ; i < N ; ++i)
        u[i] = rand() % PRIME;

    printf("--- AVANT ---\n");
    print_matrix(M);
    printf("u = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",u[i]);
    
    tic = wtime();
    res = linear_solv(M, u);
    toc = wtime();

    printf("\n--- APRES ---\nres = ");
    for(i = 0 ; i < N ; ++i)
        printf("%u ",res[i]);

    time = toc - tic;
    printf("\ntime = %f\n", time);

    free_matrix(M);
    free(u);
    free(res);   
}

void test_pluq_inversion(u32 PRIME, u32 N){
    printf("--- TEST PLUQ_INVERSION ---\n");
    Matrix *M = init_matrix(PRIME, N);
    Matrix *res;
    int i;
    double tic, toc, time;

    for (i = 0 ; i < N * N ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- M ---\n");
    print_matrix(M);

    tic = wtime();
    res = pluq_inversion(M);
    toc = wtime();

    printf("--- M^-1 ---\n");
    print_matrix(res);

    printf("--- VERIF ---\n");
    Matrix *tmp_strassen = strassen_mul(M, res);
    print_matrix(tmp_strassen);
    free_matrix(tmp_strassen);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
    free_matrix(res);
}

void test_strassen_inversion(u32 PRIME, u32 N){
    printf("--- TEST STRASSEN_INVERSION ---\n");
    Matrix *M = init_matrix(PRIME, N);
    Matrix *res;
    int i;
    double tic, toc, time;

    for (i = 0 ; i < N * N ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- M ---\n");
    print_matrix(M);

    tic = wtime();
    res = strassen_inversion(M);
    toc = wtime();

    printf("--- M^-1 ---\n");
    print_matrix(res);

    printf("--- VERIF ---\n");
    Matrix *tmp_strassen = strassen_mul(M, res);
    print_matrix(tmp_strassen);
    free_matrix(tmp_strassen);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
    free_matrix(res);
}

void test_strassen_mul_inversion(u32 PRIME, u32 N){
    printf("--- TEST STRASSEN_MUL_INVERSION ---\n");
    Matrix *M = init_matrix(PRIME, N);
    Matrix *res;
    int i;
    double tic, toc, time;

    for (i = 0 ; i < N * N ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- M ---\n");
    print_matrix(M);

    tic = wtime();
    res = strassen_mul_inversion(M);
    toc = wtime();

    printf("--- M^-1 ---\n");
    print_matrix(res);

    printf("--- VERIF ---\n");
    Matrix *tmp_strassen = strassen_mul(M, res);
    print_matrix(tmp_strassen);
    free_matrix(tmp_strassen);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
    free_matrix(res);
}

void test_strassen_optimized_inversion(u32 PRIME, u32 N){
    printf("--- TEST STRASSEN_OPTIMIZED_INVERSION ---\n");
    Matrix *M = init_matrix(PRIME, N);
    Matrix *res;
    int i;
    double tic, toc, time;

    for (i = 0 ; i < N * N ; ++i) 
        M->m[i] = rand() % PRIME;

    printf("--- M ---\n");
    print_matrix(M);

    tic = wtime();
    res = strassen_optimized_inversion(M);
    toc = wtime();

    printf("--- M^-1 ---\n");
    print_matrix(res);

    printf("--- VERIF ---\n");
    Matrix *tmp_strassen = strassen_mul(M, res);
    print_matrix(tmp_strassen);
    free_matrix(tmp_strassen);

    time = toc - tic;
    printf("time = %f\n", time);

    free_matrix(M);
    free_matrix(res);
}

void test_miller_rabin(u32 PRIME){
    printf("--- TEST MILLER_RABIN ---\n");
    u32 a = 1069639009;
    // unsigned int i;
    // for(i = 0 ; i < 10 ; i++){
    //     a = rand() % (PRIME + 1);
    //     if(miller_rabin(a)){
    //         printf("%u est un nombre premier\n", a);
    //     }
    // }

    if(miller_rabin(a)){
        printf("%u est un nombre premier\n", a);
    }
}
