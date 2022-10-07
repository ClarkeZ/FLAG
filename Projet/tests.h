/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#ifndef TESTS_H_
#define TESTS_H_

#include "algo.h"

/* ---------- TEST BASE.C ---------- */

void test_add(u32 PRIME);

void test_sub(u32 PRIME);

void test_mul(u32 PRIME);

void test_mulmod(u32 PRIME);

void test_invmod(u32 PRIME);

void test_exp_mod(u32 PRIME);

/* ---------- TEST MATRIX.C ---------- */

void test_init_matrix(u32 PRIME, u32 N);

void test_init_eye(u32 PRIME, u32 N);

void test_copy_matrix(u32 PRIME, u32 N);

void test_matrix_add(u32 PRIME, u32 N);

void test_matrix_sub(u32 PRIME, u32 N);

void test_matrix_mul(u32 PRIME, u32 N);

void test_strassen_mul(u32 PRIME, u32 N);

void test_matrix_mul_coef(u32 PRIME, u32 N);

void test_matrix_mul_optimized(u32 PRIME, u32 N);

void test_transpose(u32 PRIME, u32 N);

/* ---------- TEST ALGO.C ---------- */

void test_swap_lines(u32 PRIME, u32 N);

void test_swap_columns(u32 PRIME, u32 N);

void test_swap_pluq_pivot(u32 PRIME, u32 N);

void test_mul_line(u32 PRIME, u32 N);

void test_multiply_line(u32 PRIME, u32 N);

void test_substract_lines(u32 PRIME, u32 N);

void test_sub_lines(u32 PRIME, u32 N);

void test_zero_in_diagonal(u32 PRIME, u32 N);

void test_solv_upper_triangular_matrix(u32 PRIME, u32 N);

void test_solv_lower_triangular_matrix(u32 PRIME, u32 N);

void test_inverse_upper_triangular_matrix(u32 PRIME, u32 N);

void test_inverse_lower_triangular_matrix(u32 PRIME, u32 N);

void test_matrix_vector_product(u32 PRIME, u32 N);

void test_pluq(u32 PRIME, u32 N);

void test_linear_solv(u32 PRIME, u32 N);

void test_pluq_inversion(u32 PRIME, u32 N);

void test_strassen_inversion(u32 PRIME, u32 N);

void test_strassen_mul_inversion(u32 PRIME, u32 N);

void test_strassen_optimized_inversion(u32 PRIME, u32 N);

void test_miller_rabin(u32 PRIME);

#endif // TESTS_H