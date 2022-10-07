/*
 * Authors : Jacques Colin & Clarke Zhou
 * 
 * v1.00 (2022-05)
 * 
 * USAGE : 
 *      $ ./main -p 65537 -s 10 -i 5
 * 
 */
#include "unit_test.h"
#include "benchmark.h"

int main(int argc, char *argv[])
{
    srand(time(NULL));
    u32 prime = 1069639009;
    u32 n = 128;
    unsigned int ite = 1;

    // Affectation des valeurs : prime et ou n et ou iteration
    if(argc == 2 || argc == 3 || argc == 5 || argc == 7){
        // Affichage aide
        if( strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0) {
            printf("Commande : \n./ [-prime p] [-size s] [-iteration i]\n");
            printf("Commande : \n./ [-p pr] [-s sz] [-i it]\n");
            return EXIT_SUCCESS;
        }

        if(strcmp(argv[1], "-prime") == 0 || strcmp(argv[1], "-p") == 0){
            prime = atoi(argv[2]);
        }
        if(strcmp(argv[1], "-size") == 0 || strcmp(argv[1], "-s") == 0){
            n = atoi(argv[2]);
        }
        if(strcmp(argv[1], "-iteration") == 0 || strcmp(argv[1], "-i") == 0){
            ite = atoi(argv[2]);
        }
        
        if(argc >= 5){
            if(strcmp(argv[3], "-prime") == 0 || strcmp(argv[3], "-p") == 0){
                prime = atoi(argv[4]);
            }
            if(strcmp(argv[3], "-size") == 0 || strcmp(argv[3], "-s") == 0){
                n = atoi(argv[4]);
            }
            if(strcmp(argv[3], "-iteration") == 0 || strcmp(argv[3], "-i") == 0){
                ite = atoi(argv[4]);
            }
            
            if(argc >= 7){
                if(strcmp(argv[5], "-prime") == 0 || strcmp(argv[5], "-p") == 0){
                    prime = atoi(argv[6]);
                }
                if(strcmp(argv[5], "-size") == 0 || strcmp(argv[5], "-s") == 0){
                    n = atoi(argv[6]);
                }
                if(strcmp(argv[5], "-iteration") == 0 || strcmp(argv[5], "-i") == 0){
                    ite = atoi(argv[6]);
                }
            }
        }
    }
    else if(argc != 1 && argc != 3 && argc != 5 && argc != 7){
        printf("argc = %d\n", argc);
        fprintf(stderr, "Nombre d'argument invalide\n");
        return EXIT_FAILURE;
    }

    if(!miller_rabin(prime)){
        printf("%u n'est pas un nombre premier\n", prime);
        return EXIT_FAILURE;
    }

    printf("\n--- Initialisation ---\n");
    printf("Nombre premier : %u\n", prime);
    printf("Taille de la matrice : %u\n", n);
    printf("Nombre d'iteration(s) : %u\n\n", ite);


    /* ---------- TESTS UNITAIRES AVEC LE TEMPS D'EXECUTION ---------- */

    test_base(prime);

    test_matrix(prime, n);

    test_algo(prime, n);


    /* ---------- BENCHMARKS ---------- */

    int choix = 0;

    printf("\n--- Choix du benchmark ---\n");
    printf(" 1 - Comparaison entre la multiplication naive et la multiplication de Strassen\n");
    printf(" 2 - Comparaison entre l'inversion PLUQ et l'inversion de Strassen\n");
    printf(" 3 - Comparaison des algorithmes d'inversions de Strassen utilisant différent algorithme de multiplication\n");
    printf("        => Multiplication naive\n");
    printf("        => Multiplication de Strassen\n");
    printf("        => Multiplication naive et Strassen (choisit l'un ou l'autre selon la taille de la matrice)\n");
    printf(" 4 - Temps d'execution de l'algorithme PLUQ\n");
    printf(" 5 - Quitter\n");

    scanf("%d", &choix);

    while(choix >= 1 || choix < 5){
        switch(choix){
            case 1 : 
                /* Comparaison entre la multiplication naive et Strassen */
                benchmark_mul_naive_vs_strassen(prime, n, ite);
            break;
            case 2 : 
                /* Comparaison entre l'inversion PLUQ et Strassen */
                benchmark_inversion_pluq_vs_strassen(prime, n, ite);
            break;
            case 3 :
                /* Comparaison de l'inversion de Strassen avec la multiplication naive, strassen et optimise*/
                benchmark_strassen_inversion_StrassenMul_vs_naive_vs_optimized(prime, n, ite);
            break;
            case 4 :
                /* Temps moyen d'execution de PLUQ */
                benchmark_pluq(prime, n, ite);
            break;
            case 5 :
                return EXIT_SUCCESS;
            default :
                printf("Veuillez choisir un nombre entre 1 et 5\n");
            break;
        }
        if(choix != 5){
            printf("\n--- Choix du benchmark ---\n");
            printf(" 1 - Comparaison entre la multiplication naive et la multiplication de Strassen\n");
            printf(" 2 - Comparaison entre l'inversion PLUQ et l'inversion de Strassen\n");
            printf(" 3 - Comparaison des algorithmes d'inversions de Strassen utilisant différent algorithme de multiplication\n");
            printf("        => Multiplication naive\n");
            printf("        => Multiplication de Strassen\n");
            printf("        => Multiplication naive et Strassen (choisit l'un ou l'autre selon la taille de la matrice)\n");
            printf(" 4 - Temps d'execution de l'algorithme PLUQ\n");
            printf(" 5 - Quitter\n");

            scanf("%d", &choix);
        }
    }
    return 0;
}