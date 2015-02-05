#!/bin/bash

# Mac:
CXX=g++-4.8
# Anything else:
#CXX=g++

if [ "$1" == "all" ]
    then
        printf "Recompiling dependencies...\n"
        $CXX -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -c -o  "utils.o" "utils.cpp"
        $CXX -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -c -o  "node.o" "node.cpp"
        $CXX -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -c -o  "tree.o" "tree.cpp"
        $CXX -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -c -o  "tree_reader.o" "tree_reader.cpp"
        $CXX -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -c -o  "tree_utils.o" "tree_utils.cpp"
        $CXX -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -c -o  "superdouble.o" "superdouble.cpp"
fi

$CXX -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 -c -o "main_test_JWB.o" "main_test_JWB.cpp"
$CXX -o "pxJWB" -O3 -ffast-math -ftree-vectorize -fopenmp -std=c++11 main_test_JWB.o ./utils.o ./node.o ./tree.o ./tree_reader.o ./tree_utils.o ./superdouble.o -llapack -lblas -lpthread -lm -lnlopt_cxx -larmadillo

