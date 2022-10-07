/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#ifndef BASE_H_
#define BASE_H_

#include <stdlib.h>
#include <stdint.h>
#include "matrix.h"

typedef uint32_t u32;
typedef uint64_t u64;

/* Addition de deux entiers */
u32 add(u32 a, u32 b, u32 n);

/* Soustraction de deux entiers */
u32 sub(u32 a, u32 b, u32 n);

/* Multiplication de deux entiers */
u32 mul(u32 a, u32 b, u32 n);

/* Multiplication modulaire */
u32 mulmod(u32 a, u32 b, u32 n);

/* Inverse modulaire d'un entier a par rapport a b */
u32 invmod(u32 a, u32 b);

/* Exponentiation modulaire */
u32 exp_mod(u32 base, u32 exposant, u32 n);

#endif // BASE_H