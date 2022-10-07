/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#include "unit_test.h"

/* ---------- TEST BASE.C ---------- */
void test_base(u32 prime){

    // test_add(prime);

    // test_sub(prime);

    // test_mul(prime);

    // test_mulmod(prime);

    // test_invmod(prime);

    // test_exp_mod(prime);
}

/* ---------- TEST MATRIX.C ---------- */
void test_matrix(u32 prime, u32 n){

    // test_init_matrix(prime, n);

    // test_init_eye(prime, n);

    // test_copy_matrix(prime, n);

    // test_matrix_add(prime, n);

    // test_matrix_sub(prime, n);

    test_matrix_mul(prime, n);

    test_strassen_mul(prime, n);

    // test_matrix_mul_coef(prime, n);

    // test_matrix_mul_optimized(prime, n);

    // test_transpose(prime, n);
}

 /* ---------- TEST ALGO.C ---------- */
void test_algo(u32 prime, u32 n){

    // test_swap_lines(prime, n);

    // test_swap_columns(prime, n);

    // test_swap_pluq_pivot(prime, n);

    // test_mul_line(prime, n);

    // test_multiply_line(prime, n);

    // test_substract_lines(prime, n);

    // test_sub_lines(prime, n);

    // test_zero_in_diagonal(prime, n);

    // test_solv_upper_triangular_matrix(prime, n);

    // test_solv_lower_triangular_matrix(prime, n);

    // test_inverse_upper_triangular_matrix(prime, n);

    // test_inverse_lower_triangular_matrix(prime, n);

    // test_matrix_vector_product(prime, n);

    test_pluq(prime, n);

    // test_linear_solv(prime, n);

    // test_pluq_inversion(prime, n);

    // test_strassen_inversion(prime, n);

    // test_strassen_mul_inversion(prime, n);

    // test_strassen_optimized_inversion(prime, n);

    // test_miller_rabin(prime);
}