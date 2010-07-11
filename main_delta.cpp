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
#include "delta.h"


int main(int argc, char * argv[]){
	TreeReader tr;

	if (argc != 4){
		cout << "usage: phyx_delta l r o" << endl;
		exit(0);
	}
	Delta delta;
	vector<double> nums = delta.delta(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
	cout << nums[0] << " " << nums[1] << endl;
	return EXIT_SUCCESS;
}
