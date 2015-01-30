#!/bin/bash
if [ "$1" == "all" ]
    then
        printf "Recompiling dependencies...\n"
        g++ -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -march=native -c -o  "utils.o" "utils.cpp"
        g++ -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -march=native -c -o  "node.o" "node.cpp"
        g++ -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -march=native -c -o  "tree.o" "tree.cpp"
        g++ -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -march=native -c -o  "tree_reader.o" "tree_reader.cpp"
        g++ -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -march=native -c -o  "tree_utils.o" "tree_utils.cpp"
        g++ -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -march=native -c -o  "superdouble.o" "superdouble.cpp"
fi

g++ -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -march=native -c -o  "main_test_JWB.o" "main_test_JWB.cpp"
g++ -o "pxJWB" -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -march=native main_test_JWB.o ./utils.o ./node.o ./tree.o ./tree_reader.o ./tree_utils.o ./superdouble.o -llapack -lblas -lpthread -lm -lnlopt_cxx -larmadillo
