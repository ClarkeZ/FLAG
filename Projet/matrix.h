/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <string.h>
#include "base.h"

typedef uint32_t u32;
typedef uint64_t u64;


typedef struct Matrix
{
    u32 *m;     // matrice  M
    u32 prime;  // Nombre premier pour le corps
    u32 n;      // Taille de la matrice, pour une matrice carree
} Matrix;

typedef struct PLUQ
{
    Matrix *P;
    Matrix *L;
    Matrix *U;
    Matrix *Q;
} PLUQ;

/* Initialise une matrice de taille n*n sur Z/primeZ */
Matrix *init_matrix(u32 prime, u32 n);

/* Initialise une matrice identite de taille n*n */
Matrix *init_eye(u32 prime, u32 n);

/* Libere la memoire */
void free_matrix(Matrix *m);

/* Libere la memoire */
void free_pluq(PLUQ* pluq);

/* Copie la matrice */
Matrix *copy_matrix(Matrix *M);

/* Affiche la matrice */
void print_matrix(Matrix *M);

/* Addition de deux matrices */
Matrix *matrix_add(Matrix* A, Matrix* B);

/* Soustraction de deux matrices */
Matrix *matrix_sub(Matrix* A, Matrix* B);

/* Multiplication naive de deux matrices */
Matrix *matrix_mul(Matrix *A, Matrix *B);

/* Multiplication de deux matrices en utilisant l'algorithme de Strassen */
Matrix *strassen_mul(Matrix *A, Matrix *X);

/* Multiplication d'une matrice par un coefficient */
Matrix *matrix_mul_coef(Matrix *A, u32 c);

/* Multiplication de deux matrices en utilisant soit l'algorithme naive soit Strassen */
Matrix *matrix_mul_optimized(Matrix *A, Matrix *B);

/* Transpose la matrice */
Matrix *transpose(Matrix *A);

#endif // matrix_H