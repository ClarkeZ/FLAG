/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#include "matrix.h"

/*
Initialise une matrice carree de taille n et remplie de 0
@param prime : un nombre premier
@param n : la taille de la matrice

@return une matrice de coefficient initialise a 0 de taille n*n 
*/
Matrix *init_matrix(u32 prime, u32 n){
    Matrix *mat = malloc(sizeof(Matrix));

    if(mat == NULL){
        printf("init_matrix : Erreur allocation mémoire struct Matrix !\n");
        return NULL;
    }

    mat->n = n;
    mat->prime = prime;
    mat->m = calloc(n*n,sizeof(u32));

    if(mat->m == NULL){
        printf("init_matrix : Erreur allocation mémoire Matrix !\n");
        return NULL;
    }
    return mat;
}

/*
Initialise une matrice identite de taille n*n
@param prime : un nombre premier
@param n : la taille de la matrice

@return une matrice identite de taille n*n 
*/
Matrix *init_eye(u32 prime, u32 n) {
    unsigned int i;
    
    Matrix *mat = malloc(sizeof(Matrix));

    if(mat == NULL){
        printf("init_matrix : Erreur allocation mémoire struct Matrix !\n");
        return NULL;
    }

    mat->n = n;
    mat->prime = prime;
    mat->m = calloc(n*n,sizeof(u32));

    if(mat->m == NULL){
        printf("init_matrix : Erreur allocation mémoire Matrix !\n");
        return NULL;
    }
    for (i = 0; i < n; ++i) {
        mat->m[i*n+i] = 1;
    }
    return mat;
}

/*
Libere l'espace memoire d'une structure Matrix
*/
void free_matrix(Matrix *M){
    free(M->m);
    free(M);
}

/*
Libere l'espace memoire d'une structure PLUQ
*/
void free_pluq(PLUQ* pluq){
    free_matrix(pluq->P);
    free_matrix(pluq->L);
    free_matrix(pluq->U);
    free_matrix(pluq->Q);
    free(pluq);
}

/*
Copie une matrice dans une nouvelle matrice
@param M : une matrice a copier

@return la matrice copie
*/
Matrix *copy_matrix(Matrix *M){
    Matrix *copy = init_matrix(M->prime, M->n);

    memcpy(copy->m, M->m, sizeof(u32) * M->n * M->n);

    return copy;
}

/*
Affichage d'une matrice
@param M : une matrice a afficher
*/
void print_matrix(Matrix *M){
    unsigned int i,j;

    for(i = 0 ; i < M->n ; i++) {
        for(j = 0 ; j < M->n ; j++) {
            printf("%u ", M->m[i * M->n + j]);
        }
        printf("\n");
    }
}

/*
Addition de deux matrices A + B
@param A : une matrice 
@param B : une matrice 

@return la matrice A + B
*/
Matrix *matrix_add(Matrix *A, Matrix *B){
    unsigned int i;

    if (A->n == B->n) {
        Matrix *mat_add = init_matrix(A->prime, A->n);
        for (i = 0 ; i < A->n * A->n ; ++i) 
            mat_add->m[i] = add(A->m[i], B->m[i], A->prime);
        return mat_add;   
    }
    else {
        printf("Addition de matrice impossible, dimensions : %u != %u", A->n, B->n);
        return NULL;
    }
}

/*
Soustraction de deux matrices A - B
@param A : une matrice 
@param B : une matrice 

@return la matrice A - B
*/
Matrix *matrix_sub(Matrix *A, Matrix *B){
    unsigned int i;

    if (A->n == B->n) {
        Matrix *mat_sub = init_matrix(A->prime, A->n);
        for (i = 0 ; i < A->n * A->n ; ++i) {
            mat_sub->m[i] = sub(A->m[i], B->m[i], A->prime);
        }
        return mat_sub;   
    }
    else {
        printf("Soustraction de matrice impossible, dimensions : %u != %u", A->n, B->n);
        return NULL;
    }
}

/*
Multiplication naive de deux matrices A * B
@param A : une matrice 
@param B : une matrice 

@return la matrice A * B
*/
Matrix *matrix_mul(Matrix *A, Matrix *B) {
    unsigned int i, j, k;

    if (A->n == B->n) {
        u32 n = A->n;
        Matrix *mat_mult = init_matrix(A->prime, A->n);
        for (i = 0 ; i < n ; i++)
            for (j = 0 ; j < n ; j++)
                for (k = 0 ; k < n ; k++) {
                    u64 x = mat_mult->m[i * n + j];
                    u64 y = A->m[k + n * i];
                    u64 z = B->m[k * n + j];
                    mat_mult->m[i * n + j] = add(x, mul(y, z, A->prime), A->prime);
                }
        return mat_mult;   
    }
    else {
        printf("Multiplication de matrice impossible, dimensions : %u != %u", A->n, B->n);
        return NULL;
    }
}

/*
Multiplication de deux matrices A * B en utilisant l'algorithme de Strassen
@param A : une matrice 
@param B : une matrice 

@return la matrice A * B
*/
Matrix *strassen_mul(Matrix *A, Matrix *X){
    unsigned int i;
    if(A->n == 1 || A->n % 2 == 1 || A->n < 30)
        return matrix_mul(A, X);

    u32 n = A->n;
    u32 prime = A->prime;
    u32 n2 = A->n / 2;

    Matrix *a = init_matrix(prime, n2);
    Matrix *b = init_matrix(prime, n2);
    Matrix *c = init_matrix(prime, n2);
    Matrix *d = init_matrix(prime, n2);

    for(i = 0 ; i < n2 ; i++){
        memcpy(a->m+(i * n2), A->m+(i * n), sizeof(u32) * n2);
        memcpy(b->m+(i * n2), A->m+(i * n + n2), sizeof(u32) * n2);
        memcpy(c->m+(i * n2), A->m+(i * n + n * n2), sizeof(u32) * n2);
        memcpy(d->m+(i * n2), A->m+(i * n + n * n2 + n2), sizeof(u32) * n2);
    }

    Matrix *x = init_matrix(prime, n2);
    Matrix *y = init_matrix(prime, n2);
    Matrix *z = init_matrix(prime, n2);
    Matrix *t = init_matrix(prime, n2);

    for(i = 0 ; i < n2 ; i++){
        memcpy(x->m+(i * n2), X->m+(i * n), sizeof(u32) * n2);
        memcpy(y->m+(i * n2), X->m+(i * n + n2), sizeof(u32) * n2);
        memcpy(z->m+(i * n2), X->m+(i * n + n * n2), sizeof(u32) * n2);
        memcpy(t->m+(i * n2), X->m+(i * n + n * n2 + n2), sizeof(u32) * n2);
    }

    Matrix *tmp_q1 = matrix_add(x, z);
    Matrix *q1 = strassen_mul(a, tmp_q1);
    free_matrix(tmp_q1);

    Matrix *tmp_q2 = matrix_add(y, t);
    Matrix *q2 = strassen_mul(d, tmp_q2);
    free_matrix(tmp_q2);

    Matrix *tmp1_q3 = matrix_sub(d, a); 
    Matrix *tmp2_q3 = matrix_sub(z, y);
    Matrix *q3 = strassen_mul(tmp1_q3, tmp2_q3);
    free_matrix(tmp1_q3);
    free_matrix(tmp2_q3);

    Matrix *tmp1_q4 = matrix_sub(b, d);
    Matrix *tmp2_q4 = matrix_add(z, t);
    Matrix *q4 = strassen_mul(tmp1_q4, tmp2_q4);
    free_matrix(tmp1_q4);
    free_matrix(tmp2_q4);

    Matrix *tmp_q5 = matrix_sub(b, a);
    Matrix *q5 = strassen_mul(tmp_q5, z);
    free_matrix(tmp_q5);

    Matrix *tmp1_q6 = matrix_sub(c, a);
    Matrix *tmp2_q6 = matrix_add(x, y);
    Matrix *q6 = strassen_mul(tmp1_q6, tmp2_q6);
    free_matrix(tmp1_q6);
    free_matrix(tmp2_q6);

    Matrix *tmp_q7 = matrix_sub(c, d);
    Matrix *q7 = strassen_mul(tmp_q7, y);
    free_matrix(tmp_q7);


    Matrix *r11 = matrix_add(q1, q5);

    Matrix *tmp1_r12 = matrix_sub(q4, q5);
    Matrix *tmp2_r12 = matrix_add(q3, tmp1_r12);
    Matrix *r12 = matrix_add(q2, tmp2_r12);
    free_matrix(tmp1_r12);
    free_matrix(tmp2_r12);

    Matrix *tmp1_r21 = matrix_sub(q6, q7);
    Matrix *tmp2_r21 = matrix_add(q3, tmp1_r21);
    Matrix *r21 = matrix_add(q1, tmp2_r21);
    free_matrix(tmp1_r21);
    free_matrix(tmp2_r21);

    Matrix *r22 = matrix_add(q2, q7);

    Matrix *res = init_matrix(A->prime, A->n);

    for(i = 0; i < n/2 ; i++){
        memcpy(res->m+(i * n), r11->m+(i * n2), sizeof(u32) * n2);
        memcpy(res->m+(i * n + n2), r12->m+(i * n2), sizeof(u32) * n2);
        memcpy(res->m+(i * n + n * n2), r21->m+(i * n2), sizeof(u32) * n2);
        memcpy(res->m+(i * n + n * n2 + n2), r22->m+(i * n2), sizeof(u32) * n2);
    }

    free_matrix(a);
    free_matrix(b);
    free_matrix(c);
    free_matrix(d);
    free_matrix(x);
    free_matrix(y);
    free_matrix(z);
    free_matrix(t);

    free_matrix(q1);
    free_matrix(q2);
    free_matrix(q3);
    free_matrix(q4);
    free_matrix(q5);
    free_matrix(q6);
    free_matrix(q7);

    free_matrix(r11);
    free_matrix(r12);
    free_matrix(r21);
    free_matrix(r22);

    return res;
}   

/*
Multiplication d'une matrice par un coefficient A * c
@param A : une matrice
@param c : le coefficient 

@return la matrice A + B
*/
Matrix *matrix_mul_coef(Matrix *A, u32 c){
    unsigned int i, j;
    Matrix *res = copy_matrix(A);

    for(i = 0 ; i < A->n ; i++){
        for(j = 0 ; j < A->n ; j++){
            res->m[i*A->n + j] = mul(c, A->m[i*A->n + j],A->prime);
        }
    }
    return res;
}

/*
Multiplication de deux matrices en utilisant soit l'algorithme naive, soit l'algorithme de Strassen
@param A : une matrice
@param B : une matrice

@return la matrice A * B
*/
Matrix *matrix_mul_optimized(Matrix *A, Matrix *B){
    if(A->n < 70)
        return matrix_mul(A, B);
    else
        return strassen_mul(A, B);
}

/*
Transpose une matrice A
@param A : une matrice

@return la matrice A transposee A^T
*/
Matrix *transpose(Matrix *A){
    unsigned int i, j;

    Matrix *transp = init_matrix(A->prime, A->n);
    for(i = 0 ; i < A->n ; i++)
        for(j = 0 ; j < A->n ; j++)
            transp->m[i * A->n + j] = A->m[j * A->n + i];

    return transp;
}