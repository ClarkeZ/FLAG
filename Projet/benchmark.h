/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include "algo.h"

/* Comparaison entre la multiplication naive et Strassen */
void benchmark_mul_naive_vs_strassen(u32 PRIME, u32 N, unsigned int ITE);

/* Comparaison entre l'inversion PLUQ et Strassen */
void benchmark_inversion_pluq_vs_strassen(u32 PRIME, u32 N, unsigned int ITE);

/* Comparaison de l'inversion de Strassen avec la multiplication naive, strassen et optimise*/
void benchmark_strassen_inversion_StrassenMul_vs_naive_vs_optimized(u32 PRIME, u32 N, unsigned int ITE);

/* Temps d'execution de PLUQ */
void benchmark_pluq(u32 PRIME, u32 N, unsigned int ITE);


#endif // BENCHMARK_H