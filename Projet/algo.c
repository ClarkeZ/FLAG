/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 * USAGE : 
 *      $ ./main -p 65537 -s 10 -i 5
 * 
 */
#include "algo.h"

double wtime(){
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

/*
Interchange deux lignes d'une matrice
@param M : une matrice
@param l1 : un numero de ligne
@param l2 : un numero de ligne
*/
void swap_lines(Matrix *M, u32 l1, u32 l2){
    if(l1 == l2)
        return;

    u32 *tmp = (u32 *) malloc(sizeof(u32) * M->n);
    memcpy(tmp,              M->m + l1 * M->n,  M->n * sizeof(u32));
    memcpy(M->m + l1 * M->n, M->m + l2 * M->n,  M->n * sizeof(u32));
    memcpy(M->m + l2 * M->n, tmp,               M->n * sizeof(u32));
    free(tmp);
}

/*
Interchange deux colonnes d'une matrice
@param M : une matrice
@param c1 : un numero de colonne
@param c2 : un numero de colonne
*/
void swap_columns(Matrix *M, u32 c1, u32 c2){
    if(c1 == c2)
        return;

    u32 tmp = 0;
    for(unsigned i = 0 ; i < M->n * M->n ; i+= M->n){
        tmp = M->m[c1 + i];
        M->m[c1 + i] = M->m[c2 + i];
        M->m[c2 + i] = tmp;
    }
}

/*
Cherche le coefficient le plus eleve et le place a la position du pivot (i.e dans la diagonale)
@param pluq : une structure pluq 
@param min : l'emplacement du pivot (min = ligne = colonne)
*/
void swap_pluq_pivot(PLUQ *pluq, unsigned int min) {
    u32 i, j, l, c, max, n;
    n = pluq->U->n;
    max = 0;
    l = min;
    c = min;

    /* Recherche le plus grand coefficient */
    for (i = min; i < n; ++i)
        for (j = min; j < n; ++j)
            if (pluq->U->m[i*n+j] > max){
                max = pluq->U->m[i*n+j];
                l = i;
                c = j;
            }

    /* SWAP */
    if(min != l){
        swap_columns(pluq->P, min, l);
        swap_lines(pluq->L, min, l);
        swap_columns(pluq->L, min, l);
        swap_lines(pluq->U, min, l);
    }
    if(min != c){
        swap_columns(pluq->U, min, c);
        swap_lines(pluq->Q, min, c);
        
    }
}

/*
Multiplication d'une ligne d'une matrice par un coefficient
@param A : une matrice
@param coef : un coefficient
@param ligne : un numero de ligne de la matrice A
*/
void mul_line(Matrix *A, u32 coef, unsigned int ligne){
    unsigned int k;

    for(k = A->n * ligne ; k < A->n * (ligne + 1) ; k++){
        A->m[k] = mul(A->m[k], coef, A->prime);
    }
}

/*
Multiplication d'une ligne d'une matrice (un vecteur) par un coefficient
@param line : une ligne de la matrice, ou un vecteur
@param coef : un coefficient
@param n : la taille de la matrice (de la ligne)
@param prime : un nombre premier

@return le vecteur multiplier par le coefficient line * coef 
*/
u32 *multiply_line(u32 *line, u32 coef, u32 n, u32 prime){
    unsigned int i;
    u32 *res = (u32 *) malloc(sizeof(u32) * n);

    for(i = 0 ; i < n ; i++)
        res[i] = mul(line[i], coef, prime);

    return res;
}

/*
Soustraction d'un vecteur par un autre vecteur
@param line1 : une ligne de la matrice, ou un vecteur
@param line2 : une ligne de la matrice, ou un vecteur
@param n : la taille de la matrice (de la ligne)
@param prime : un nombre premier

@return le vecteur line soustrait par le vecteur line2 line1 - line2
*/
u32 *substract_lines(u32 *line1, u32 *line2, u32 n, u32 prime){
    unsigned int i;
    u32 *res = (u32 *) malloc(sizeof(u32) * n);

    for(i = 0 ; i < n ; i++)
        res[i] = sub(line1[i], line2[i], prime);

    return res;
}

/*
Soustration d'une ligne i par une ligne j de la matrice
@param A : une matrice
@param i : un numero de ligne de la matrice A
@param j : un numero de ligne de la matrice A
*/
void sub_lines(Matrix *A, unsigned int i, unsigned int j){
    unsigned int k;
    i = i * A->n;
    j = j * A->n;

    for(k = 0 ; k < A->n ; k++) {
        A->m[i + k] = sub(A->m[i + k], A->m[j + k], A->prime);
    }
}

/*
Verifie si la diagonale de la matrice possede un coefficient nul
@param A : une matrice

@return true s'il existe un coefficient nul, false sinon
*/
bool zero_in_diagonal(Matrix *A){
    unsigned int i;

    for(i = 0 ; i < A->n ; i++)
        if(A->m[i * A->n + i] == 0)
            return true;   
    
    return false;
}

/*
Resolution lineaire d'une matrice triangulaire superieure Ux = v
@paraim U : une matrice superieure
@param v : un vecteur, le resultat

@return un vecteur x resolvant Ux = v
*/
u32 *solv_upper_triangular_matrix(Matrix *U, u32 *v){
    int i, j;
    
    u32 *x = (u32 *) malloc(sizeof(u32) * U->n);
    memcpy(x, v, sizeof(u32) * U->n);

    for(i = U->n - 1 ; i >= 0 ; i--){
        for(j = U->n - 1 ; j > i ; j--)
            x[i] = (sub(x[i], mul(x[j], U->m[i * U->n + j], U->prime), U->prime)) ;
    
        x[i] = mul(invmod(U->m[i * U->n + i], U->prime), x[i], U->prime);
    }
    return x;
}

/*
Resolution lineaire d'une matrice triangulaire inferieure Lx = v
@paraim L : une matrice superieure
@param v : un vecteur, le resultat

@return un vecteur x resolvant Lx = v
*/
u32 *solv_lower_triangular_matrix(Matrix *L, u32 *v){
    unsigned int i, j;
    
    u32 *x = (u32 *) malloc(sizeof(u32) * L->n);
    memcpy(x, v, sizeof(u32) * L->n);

    for(i = 0 ; i < L->n ; i++){
        for(j = 0 ; j < i ; j++)
            x[i] = (sub(x[i], mul(x[j], L->m[i * L->n + j], L->prime), L->prime)) ;

        x[i] = mul(invmod(L->m[i * L->n + i], L->prime), x[i], L->prime);
    }
    return x;
}

/*
Inversion d'une matrice triangulaire superieure
@param U : une matrice

@return une matrice inverse U^-1
*/
Matrix *inverse_upper_triangular_matrix(Matrix *U){
    if(zero_in_diagonal(U)){
        fprintf(stderr,"#inverse_upper_triangular_matrix : Inversion de la matrice U impossible car det(A) = 0\n");
        // print_matrix(U);
        free_matrix(U);
        exit(0);
    }
    int i, j;
    u32 coef = 0;
    u32 n = U->n;
    u32 prime = U->prime;
    u32 *tmp1_mat, *tmp2_mat, *tmp3_mat;
    u32 *tmp1_id, *tmp2_id, *tmp3_id;
    
    Matrix *id = init_eye(prime, n);
    Matrix *U_copy = copy_matrix(U);

    for(i = U->n - 1 ; i >= 0 ; i--){
        coef = invmod(U_copy->m[i * n + i], prime);

        tmp1_mat = multiply_line(U_copy->m+(i * n), coef, n, prime);
        memcpy(U_copy->m+(i*n), tmp1_mat, sizeof(u32)*n);

        tmp1_id = multiply_line(id->m+(i * n), coef, n, prime);
        memcpy(id->m+(i*n), tmp1_id, sizeof(u32)*n);

        for(j = 0 ; j < i ; j++){
            tmp2_mat = multiply_line(tmp1_mat, U_copy->m[j*n+i], n, prime);
            tmp2_id = multiply_line(tmp1_id, U_copy->m[j*n+i], n, prime);

            tmp3_mat = substract_lines(U_copy->m+(j*n), tmp2_mat, n, prime);
            tmp3_id = substract_lines(id->m+(j*n), tmp2_id, n, prime);

            memcpy(U_copy->m+(j*n), tmp3_mat, sizeof(u32)*n);
            memcpy(id->m+(j*n), tmp3_id, sizeof(u32)*n);
            
            free(tmp2_mat);
            free(tmp3_mat);
            free(tmp2_id);
            free(tmp3_id);
        }
        free(tmp1_mat);
        free(tmp1_id);
    }
    free_matrix(U_copy);

    return id;
}

/*
Inversion d'une matrice triangulaire inferieure
@param L : une matrice

@return une matrice inverse L^-1
*/
Matrix *inverse_lower_triangular_matrix(Matrix *L){
    u32 i, j;
    u32 coef = 0;
    u32 n = L->n;
    u32 prime = L->prime;
    u32 *tmp1_mat, *tmp2_mat, *tmp3_mat; 
    u32 *tmp1_id, *tmp2_id, *tmp3_id;
    
    Matrix *id = init_eye(prime, n);
    Matrix *L_copy = copy_matrix(L);

    for(i = 0 ; i < n ; i++){
        coef = invmod(L_copy->m[i * n + i], prime);

        tmp1_mat = multiply_line(L_copy->m+(i * n), coef, n, prime);
        memcpy(L_copy->m+(i*n), tmp1_mat, sizeof(u32)*n);

        tmp1_id = multiply_line(id->m+(i * n), coef, n, prime);
        memcpy(id->m+(i*n), tmp1_id, sizeof(u32)*n);

        for (j = i + 1 ; j < n ; j++) {
            tmp2_mat = multiply_line(tmp1_mat, L_copy->m[j*n+i], n, prime);
            tmp2_id = multiply_line(tmp1_id, L_copy->m[j*n+i], n, prime);

            tmp3_mat = substract_lines(L_copy->m+(j*n), tmp2_mat, n, prime);
            tmp3_id = substract_lines(id->m+(j*n), tmp2_id, n, prime);

            memcpy(L_copy->m+(j*n), tmp3_mat, sizeof(u32)*n);
            memcpy(id->m+(j*n), tmp3_id, sizeof(u32)*n);

            free(tmp2_mat);
            free(tmp3_mat);
            free(tmp2_id);
            free(tmp3_id);
        }
        free(tmp1_mat);
        free(tmp1_id);
    }
    free_matrix(L_copy);

    return id;
}

/*
Produit d'une matrice avec un vecteur
@param A : une matrice
@param v : un vecteur

@return le produit matrice-vecteur A * v
*/
u32 *matrix_vector_product(Matrix *A, u32 *v){
    u32 i, j;
    u32 *res = calloc(A->n, sizeof(u32));

    for(i = 0 ; i < A->n ; i++)
        for(j = 0 ; j < A->n ; j++)
            res[i] = add(res[i], mul(A->m[i * A->n + j], v[j], A->prime), A->prime);

    return res;
}

/*
Decomposition PLUQ d'une matrice A
@param A : une matrice

@return la decomposition PLUQ de la matrice A
*/
PLUQ *pluq(Matrix *A) {
    /* Hypothèse : Les coefs de la matrice A sont mod prime */
    unsigned int i,j;
    u32 prime = A->prime;
    u32 n = A->n;
    u32 coef = 0;
    u32 *tmp1, *tmp2, *tmp3;
    
    PLUQ *pluq = (PLUQ*) malloc(sizeof(PLUQ));
    pluq->P = init_eye(prime, n);
    pluq->L = init_eye(prime, n);
    pluq->U = init_matrix(prime, n);
    pluq->Q = init_eye(prime, n);
    
    memcpy(pluq->U->m, A->m, sizeof(u32) * n * n); // Use copy_matrix(Matrix)

    for(i = 0 ; i < n ; ++i){
        swap_pluq_pivot(pluq, i);
        coef = invmod(pluq->U->m[i * n + i], prime);
        tmp1 = multiply_line(pluq->U->m+(i * n), coef, n, prime);

        for(j = i + 1 ; j < n ; j++){
            pluq->L->m[j*n+i] = add(pluq->L->m[j*n+i], mul(coef,pluq->U->m[j*n+i],prime), prime);
            /* Echelonnage */
            tmp2 = multiply_line(tmp1, pluq->U->m[j*n+i], n, prime);
            tmp3 = substract_lines(pluq->U->m+(j*n), tmp2, n, prime);
            memcpy(pluq->U->m+(j*n), tmp3, sizeof(u32)*n);

            free(tmp2);
            free(tmp3);
        }
        free(tmp1);
    }    
    return pluq;
}

/*
Resolution d'un systeme lineaire Ax = b utilisant la decomposition PLUQ
@param A : une matrice
@param b : un vecteur

@return un vecteur x qui est solution du systeme lineaire Ax = b
*/
u32 *linear_solv(Matrix *A, u32 *b){
    PLUQ *pluq_decomp = pluq(A);

    if(zero_in_diagonal(pluq_decomp->U)){
        fprintf(stderr,"#linear_solv : Ax = b impossible car det(A) = 0\n");
        // print_matrix(A);
        free_pluq(pluq_decomp);
        free_matrix(A);
        free(b);
        exit(0);
    }
    Matrix *transp_p = transpose(pluq_decomp->P);
    
    u32 *ptb = matrix_vector_product(transp_p, b);

    u32 *y = solv_lower_triangular_matrix(pluq_decomp->L, ptb);

    u32 *x_prime = solv_upper_triangular_matrix(pluq_decomp->U, y);

    Matrix *transp_q = transpose(pluq_decomp->Q);

    u32 *res = matrix_vector_product(transp_q, x_prime);

    free_pluq(pluq_decomp);
    free_matrix(transp_p);
    free_matrix(transp_q);
    free(ptb);
    free(y);
    free(x_prime);

    return res;     
}

/*
Inversion d'une matrice A en utilisant la decomposition PLUQ
@param A : une matrice

@return une matrice inverse A^-1
*/
Matrix *pluq_inversion(Matrix *A){
    PLUQ *pluq_res = pluq(A);
    if(zero_in_diagonal(pluq_res->U)){
        // printf("--- pluq inversion A ---");
        // print_matrix(A);
        fprintf(stderr,"#pluq_inversion : Inversion de la matrice impossible car det = 0\n");
        free_pluq(pluq_res);
        free_matrix(A);
        exit(0);
    }
    PLUQ *tmp = (PLUQ*) malloc(sizeof(PLUQ));
    tmp->P = transpose(pluq_res->P);
    tmp->Q = transpose(pluq_res->Q);
    tmp->L = inverse_lower_triangular_matrix(pluq_res->L);
    tmp->U = inverse_upper_triangular_matrix(pluq_res->U);
    Matrix *tmp1 = matrix_mul(tmp->Q, tmp->U);
    Matrix *tmp2 = matrix_mul(tmp1, tmp->L);
    Matrix *inv_pluq = matrix_mul(tmp2, tmp->P); 
    
    free_matrix(tmp1);
    free_matrix(tmp2);
    free_pluq(tmp);
    free_pluq(pluq_res);
    return inv_pluq;
}

/*
Inversion d'une matrice A en utilisant l'algorithme de Strassen avec la multiplication naive
@param A : une matrice

@return une matrice inverse A^-1
*/
Matrix *strassen_inversion(Matrix *M){
    unsigned int i;
    if(M->n == 1){
        Matrix *Res = init_matrix(M->prime, 1);
        Res->m[0] = invmod(M->m[0], M->prime);
        return Res;
    }
    if(M->n & 1) // (M->n & 1) = M->n % 2 == 1 
        return pluq_inversion(M);
    
    u32 n = M->n;
    u32 n2 = n / 2; 

    u32 prime = M->prime;
    
    Matrix *A = init_matrix(prime, n2);
    Matrix *B = init_matrix(prime, n2);
    Matrix *C = init_matrix(prime, n2);
    Matrix *D = init_matrix(prime, n2);

    for(i = 0 ; i < n2 ; i++){
        memcpy(A->m+(i * n2), M->m+(i * n), sizeof(u32) * n2);
        memcpy(B->m+(i * n2), M->m+(i * n + n2), sizeof(u32) * n2);
        memcpy(C->m+(i * n2), M->m+(i * n + n * n2), sizeof(u32) * n2);
        memcpy(D->m+(i * n2), M->m+(i * n + n * n2 + n2), sizeof(u32) * n2);
    }
    Matrix *E;
    if((n2 & 1) && (n2 != 1)) // (n2 & 1) = n2 % 2 == 1 
        E = pluq_inversion(A);
    else 
        E = strassen_inversion(A);

    Matrix *tmp_Z1 = matrix_mul(E, B);
    Matrix *tmp_Z2 = matrix_mul(C, tmp_Z1);
    Matrix *Z = matrix_sub(D, tmp_Z2);

    free_matrix(tmp_Z1);
    free_matrix(tmp_Z2);
    
    Matrix *T;
    if((n2 & 1) && (n2 != 1))
        T = pluq_inversion(Z);
    else
        T = strassen_inversion(Z);
    

    Matrix *tmp_invA1 = matrix_mul(C, E);
    Matrix *tmp_invA2 = matrix_mul(T, tmp_invA1);
    Matrix *tmp_invA3 = matrix_mul(B, tmp_invA2);
    Matrix *tmp_invA4 = matrix_mul(E, tmp_invA3);
    Matrix *inv_A = matrix_add(E, tmp_invA4);

    free_matrix(tmp_invA1);
    free_matrix(tmp_invA2);
    free_matrix(tmp_invA3);
    free_matrix(tmp_invA4);

    Matrix *mat_0 = init_matrix(prime, n2);
    Matrix *tmp_invB1 = matrix_mul(B, T);
    Matrix *tmp_invB2 = matrix_mul(E, tmp_invB1);
    Matrix *inv_B = matrix_sub(mat_0, tmp_invB2);

    free_matrix(tmp_invB1);
    free_matrix(tmp_invB2);

    Matrix *tmp_invC1 = matrix_mul(C, E);
    Matrix *tmp_invC2 = matrix_mul(T, tmp_invC1);
    Matrix *inv_C = matrix_sub(mat_0, tmp_invC2);

    free_matrix(tmp_invC1);
    free_matrix(tmp_invC2);
    free_matrix(mat_0);
    
    Matrix *inv = init_matrix(prime, n);

    for(i = 0; i < n/2 ; i++){
        memcpy(inv->m+(i * n), inv_A->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n2), inv_B->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n * n2), inv_C->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n * n2 + n2), T->m+(i * n2), sizeof(u32) * n2);
    }
    free_matrix(A);
    free_matrix(B);
    free_matrix(C);
    free_matrix(D);

    free_matrix(E);
    free_matrix(Z);
    free_matrix(T);

    free_matrix(inv_A);
    free_matrix(inv_B);
    free_matrix(inv_C);

    return inv;
}

/*
Inversion d'une matrice A en utilisant l'algorithme et la multiplication de Strassen
@param A : une matrice

@return une matrice inverse A^-1
*/
Matrix *strassen_mul_inversion(Matrix *M){
    unsigned int i;
    if(M->n == 1){
        Matrix *Res = init_matrix(M->prime, 1);
        Res->m[0] = invmod(M->m[0], M->prime);
        return Res;
    }
    if(M->n & 1)
        return pluq_inversion(M);

    u32 n = M->n;
    u32 n2 = n / 2; 

    u32 prime = M->prime;
    
    Matrix *A = init_matrix(prime, n2);
    Matrix *B = init_matrix(prime, n2);
    Matrix *C = init_matrix(prime, n2);
    Matrix *D = init_matrix(prime, n2);

    for(i = 0 ; i < n/2 ; i++){
        memcpy(A->m+(i * n2), M->m+(i * n), sizeof(u32) * n2);
        memcpy(B->m+(i * n2), M->m+(i * n + n2), sizeof(u32) * n2);
        memcpy(C->m+(i * n2), M->m+(i * n + n * n2), sizeof(u32) * n2);
        memcpy(D->m+(i * n2), M->m+(i * n + n * n2 + n2), sizeof(u32) * n2);
    }
    Matrix *E;
    if((n2 & 1) && (n2 != 1))
        E = pluq_inversion(A);
    else 
        E = strassen_mul_inversion(A);
    Matrix *tmp_Z1 = strassen_mul(E, B);
    Matrix *tmp_Z2 = strassen_mul(C, tmp_Z1);

    Matrix *Z = matrix_sub(D, tmp_Z2);

    free_matrix(tmp_Z1);
    free_matrix(tmp_Z2);
    
    Matrix *T;
    if((n2 & 1) && (n2 != 1))
        T = pluq_inversion(Z);
    else
        T = strassen_mul_inversion(Z);
    
    Matrix *inv = init_matrix(prime, n);
    Matrix *tmp_invA1 = strassen_mul(C, E);
    Matrix *tmp_invA2 = strassen_mul(T, tmp_invA1);
    Matrix *tmp_invA3 = strassen_mul(B, tmp_invA2);
    Matrix *tmp_invA4 = strassen_mul(E, tmp_invA3);
    Matrix *inv_A = matrix_add(E, tmp_invA4);

    free_matrix(tmp_invA1);
    free_matrix(tmp_invA2);
    free_matrix(tmp_invA3);
    free_matrix(tmp_invA4);

    Matrix *mat_0 = init_matrix(prime, n2);
    Matrix *tmp_invB1 = strassen_mul(B, T);
    Matrix *tmp_invB2 = strassen_mul(E, tmp_invB1);
    Matrix *inv_B = matrix_sub(mat_0, tmp_invB2);

    free_matrix(tmp_invB1);
    free_matrix(tmp_invB2);

    Matrix *tmp_invC1 = strassen_mul(C, E);
    Matrix *tmp_invC2 = strassen_mul(T, tmp_invC1);
    Matrix *inv_C = matrix_sub(mat_0, tmp_invC2);

    free_matrix(tmp_invC1);
    free_matrix(tmp_invC2);
    free_matrix(mat_0);

    for(i = 0; i < n/2 ; i++){
        memcpy(inv->m+(i * n), inv_A->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n2), inv_B->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n * n2), inv_C->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n * n2 + n2), T->m+(i * n2), sizeof(u32) * n2);
    }
    free_matrix(A);
    free_matrix(B);
    free_matrix(C);
    free_matrix(D);

    free_matrix(E);
    free_matrix(Z);
    free_matrix(T);

    free_matrix(inv_A);
    free_matrix(inv_B);
    free_matrix(inv_C);

    return inv;
}

/*
Inversion d'une matrice A en utilisant l'algorithme de Strassen avec soit la multiplication naive 
soit la multiplication de Strassen
@param A : une matrice

@return une matrice inverse A^-1
*/
Matrix *strassen_optimized_inversion(Matrix *M){
    unsigned int i;
    if(M->n == 1){
        Matrix *Res = init_matrix(M->prime, 1);
        Res->m[0] = invmod(M->m[0], M->prime);
        return Res;
    }
    if((M->n & 1)) // ~(M->n & 1) = M->n % 2 == 1 
        return pluq_inversion(M);

    u32 n = M->n;
    u32 n2 = n / 2; // On suppose que n est pair

    u32 prime = M->prime;
    
    Matrix *A = init_matrix(prime, n2);
    Matrix *B = init_matrix(prime, n2);
    Matrix *C = init_matrix(prime, n2);
    Matrix *D = init_matrix(prime, n2);

    for(i = 0 ; i < n/2 ; i++){
        memcpy(A->m+(i * n2), M->m+(i * n), sizeof(u32) * n2);
        memcpy(B->m+(i * n2), M->m+(i * n + n2), sizeof(u32) * n2);
        memcpy(C->m+(i * n2), M->m+(i * n + n * n2), sizeof(u32) * n2);
        memcpy(D->m+(i * n2), M->m+(i * n + n * n2 + n2), sizeof(u32) * n2);
    }
    Matrix *E;
    if((n2 & 1) && (n2 != 1))
        E = pluq_inversion(A);
    else 
        E = strassen_optimized_inversion(A);

    Matrix *tmp_Z1 = matrix_mul_optimized(E, B);
    Matrix *tmp_Z2 = matrix_mul_optimized(C, tmp_Z1);
    Matrix *Z = matrix_sub(D, tmp_Z2);

    free_matrix(tmp_Z1);
    free_matrix(tmp_Z2);
    
    Matrix *T;
    if((n2 & 1) && (n2 != 1))
        T = pluq_inversion(Z);
    else
        T = strassen_optimized_inversion(Z);
    
    Matrix *inv = init_matrix(prime, n);

    Matrix *tmp_invA1 = matrix_mul_optimized(C, E);
    Matrix *tmp_invA2 = matrix_mul_optimized(T, tmp_invA1);
    Matrix *tmp_invA3 = matrix_mul_optimized(B, tmp_invA2);
    Matrix *tmp_invA4 = matrix_mul_optimized(E, tmp_invA3);
    Matrix *inv_A = matrix_add(E, tmp_invA4);

    free_matrix(tmp_invA1);
    free_matrix(tmp_invA2);
    free_matrix(tmp_invA3);
    free_matrix(tmp_invA4);

    Matrix *mat_0 = init_matrix(prime, n2);
    Matrix *tmp_invB1 = matrix_mul_optimized(B, T);
    Matrix *tmp_invB2 = matrix_mul_optimized(E, tmp_invB1);
    Matrix *inv_B = matrix_sub(mat_0, tmp_invB2);

    free_matrix(tmp_invB1);
    free_matrix(tmp_invB2);

    Matrix *tmp_invC1 = matrix_mul_optimized(C, E);
    Matrix *tmp_invC2 = matrix_mul_optimized(T, tmp_invC1);
    Matrix *inv_C = matrix_sub(mat_0, tmp_invC2);

    free_matrix(tmp_invC1);
    free_matrix(tmp_invC2);
    free_matrix(mat_0);
    
    for(i = 0; i < n/2 ; i++){
        memcpy(inv->m+(i * n), inv_A->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n2), inv_B->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n * n2), inv_C->m+(i * n2), sizeof(u32) * n2);
        memcpy(inv->m+(i * n + n * n2 + n2), T->m+(i * n2), sizeof(u32) * n2);
    }
    free_matrix(A);
    free_matrix(B);
    free_matrix(C);
    free_matrix(D);
    free_matrix(E);
    free_matrix(Z);
    free_matrix(T);

    free_matrix(inv_A);
    free_matrix(inv_B);
    free_matrix(inv_C);
    return inv;
}

/*
Test de primalite de Miller Rabin (Algo probabiliste)
@param n : un entier

@return Vrai si n est probablement premier, faux sinon
*/
bool miller_rabin(u32 n){
    int i;
    int k = 15;
    u32 s; 

    if (n < 2)
        return 0;

    if (n != 2 && !(n & 1))
        return false;

    s = n - 1;
    while(!(s & 1)){
        s /= 2;
    }

    for (i = 0; i < k; i++)
    {
        u32 a = rand() % (n - 1) + 1, temp = s;
        u32 mod = exp_mod(a, temp, n);
        while (temp != n - 1 && mod != 1 && mod != n - 1)
        {
            mod = mulmod(mod, mod, n);
            temp *= 2;
        }
        if (mod != n - 1 && !(temp & 1))
        {
            return false;
        }
    }
    return true;
}