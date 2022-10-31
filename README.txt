Pour compiler le programme il suffit de lancer la commande suivante : 
$ make

puis lancer l’exécutable en faisant la commande:
./main [-prime p] [-size s] [-iteration i]
ou
./main [-p pr] [-s sz] [-i it]


Si aucun argument n’est mis, ou si des arguments sont manquants alors les valeurs par défaut seront :
	- prime : 1069639009
	- size : 128
	- iteration : 1

Exemple:
./main -p 65537 -s 100 -i 20
./main -p 11 -s 15
./main -s 20

Si vous avez un doute, vous pouvez utiliser la commande :

$ ./main -h
ou
$ ./main -help

Une fois exécuté, le programme vous demandera de choisir quel benchmark lancer 




Pour exécuter les tests unitaires il faut décommenter des fonctions dans le fichier unit_test.c et exécuter le programme.

