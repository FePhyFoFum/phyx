/*
 * main_mrca.cpp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "bd_sim.h"


int main(int argc, char * argv[]){
	TreeReader tr;

	if (argc != 5){
		cout << "usage: phyx_bd extant time birth death" << endl;
		exit(0);
	}
	BirthDeathSimulator bd(atof(argv[1]),atof(argv[2]),atof(argv[3]),atof(argv[4]));
	Tree * bdtr = bd.make_tree(false);
	cout << bdtr->getRoot()->getNewick(true)<<";"<<endl;
	return EXIT_SUCCESS;
}
