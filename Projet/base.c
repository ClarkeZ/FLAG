/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 */
#include "base.h"

/*
Addition de deux entiers modulo n
@param a : un entier
@param b : un entier
@param n : un modulo

@return un entier (a + b) % n
*/
u32 add(u32 a, u32 b, u32 n) {
    /* Tester temps des casts */
    if (a + b < n) {
        return a+b;
    }
    u64 a64 = a;
    u64 b64 = b;
    return (a64 + b64) % n;
}

/*
Soustraction de deux entiers modulo n
@param a : un entier
@param b : un entier
@param n : un modulo

@return un entier (a - b) % n
*/
u32 sub(u32 a, u32 b, u32 n){
    if (a >= b){
        return (a - b) % n;
    }
    else {
        // return ((n + a%n) - b%n)%n;
        return n - ((b - a) % n );
    }
}

/*
Multiplication de deux entiers modulo n
@param a : un entier
@param b : un entier
@param n : un modulo

@return un entier (a * b) % n
*/
u32 mul(u32 a, u32 b, u32 n){
    if(n <= 65537){
        return (a * b) % n;
    }
    else{
        u64 a64 = a;
        u64 b64 = b;
        return (a64 * b64) % n;
    }
}

/*
Multiplication modulaire et prend en compte le depassement (overflow)
@param a : un entier
@param b : un entier
@param n : un entier

@return un entier (a * b) mod n
*/
u32 mulmod(u32 a, u32 b, u32 n)
{
    u32 x = 0;
    u32 y = a % n;
    while (b > 0)
    {
        if (b & 1)
        {    
            x = (x + y) % n;
        }
        y = (y * 2) % n;
        b >>= 1;
    }
    return x % n;
}

/*
Inverse modulaire d'un entier
@param a : un entier
@param b : un entier

@return l'inverse modulaire a^(-1)[b]
*/
u32 invmod(u32 a, u32 b) {
    long int t = 0;  
    long int nt = 1;  
    long int r = b;  
    long int nr = a % b;
    while (nr != 0) {
        long int q = r / nr;
        long int tmp = nt;  
        nt = t - q*nt;  
        t = tmp;
        tmp = nr;  
        nr = r - q*nr;  
        r = tmp;
    }
    if (t < 0)
        t += b;
    return (u32) t;
}

/*
Exponentiation modulaire
@param base : un entier, la base
@param exposant : un entier, la puissance auquel on veut elever la base
@param n : un entier, le modulo

@return un entier base^exposant mod n
*/
u32 exp_mod(u32 base, u32 exposant, u32 n)
{
    u64 x = 1;
    u64 y = base;

    while (exposant > 0)
    {
        if (exposant & 1)
            x = (x * y) % n;
        y = (y * y) % n;
        exposant >>= 1;
    }
    return x % n;
}