/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#ifndef ALGO_H_
#define ALGO_H_

#include <stdbool.h>
#include <sys/time.h>
#include <time.h>

#include "base.h"

double wtime();

/* Interchange deux lignes d'une matrice */
void swap_lines(Matrix *M, u32 l1, u32 l2);

/* Interchange deux colonnes d'une matrice */
void swap_columns(Matrix *M, u32 c1, u32 c2);

/* Recupere le plus grand coefficient et le place a la position du pivot */
void swap_pluq_pivot(PLUQ *pluq, unsigned int min);

/* Multiplication d'une ligne d'une matrice par un coefficient */
void mul_line(Matrix *A, u32 coef, unsigned int i);

/* Multiplication d'une ligne d'une matrice par un coefficient */
u32 *multiply_line(u32 *line, u32 coef, u32 n, u32 prime);

/* Soustraction d'une ligne de matrice par une autre ligne */
u32 *substract_lines(u32 *line1, u32 *line2, u32 n, u32 prime);

/* Soustraction d'une ligne (i) par une autre ligne (j) de la Matrice */
void sub_lines(Matrix *A, unsigned int i, unsigned int j);

/* Verifie si la diagonale possede un coefficient nul */
bool zero_in_diagonal(Matrix *A);

/* Resolution de la matrice triangulaire superieure */
u32 *solv_upper_triangular_matrix(Matrix *U, u32 *v);

/* Resolution de la matrice triangulaire inferieure */
u32 *solv_lower_triangular_matrix(Matrix *L, u32 *v);

/* Inversion d'une matrice triangulaire superieure */
Matrix *inverse_upper_triangular_matrix(Matrix *U);

/* Inversion d'une matrce triangulaire inferieure */
Matrix *inverse_lower_triangular_matrix(Matrix *L);

/* Produit matrice vecteur */
u32 *matrix_vector_product(Matrix *A, u32 *v);

/* decomposition PLUQ */
PLUQ* pluq(Matrix *A);

/* Resolution d'un systeme lineaire en utilisant la decomposition PLUQ */
u32 *linear_solv(Matrix *A, u32 *b);

/* Inversion d'une matrice en utilisant PLUQ */
Matrix *pluq_inversion(Matrix *A);

/* Inversion d'une matrice en utilisant l'algorithme de Strassen */
Matrix *strassen_inversion(Matrix *M);

/* Inversion d'une matrice en utilisant l'algorithme et la multiplication de Strassen */
Matrix *strassen_mul_inversion(Matrix *M);

/* Inversion d'une matrice en utilisant l'algorithme de Strassen et la multiplication naive ou Strassen */
Matrix *strassen_optimized_inversion(Matrix *M);

/* Test de primalite de Miller Rabin */
bool miller_rabin(u32 n);

#endif // ALGO_H