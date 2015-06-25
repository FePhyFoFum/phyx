/*
 * main_seqgen.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: joe
 */


#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <map>
#include <iterator>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "seqgen.h"
#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "node.h"
#include "tree.h"
#include "tree_reader.h"

//g++ -std=c++11 branch_segment.cpp seqgen.cpp node.cpp tree.cpp tree_reader.cpp main_seqgen.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp -o test
//g++ -std=c++11 branch_segment.cpp seqgen.cpp node.cpp tree.cpp tree_reader.cpp main_seqgen.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp -larmadillo  -llapack -lblas -lpthread -lm -o test

int main(){


	TreeReader TreeUtil;
	Tree * tree;
	SEQGEN SGen;
	double length;
	string path, newick, line;
	int sequence_length = 0;
	ifstream readline;
	path = ("TestFiles/RAxML_bestTree.drosML");
	readline.open(path.c_str());
	if (readline.is_open()){
		while (getline(readline, line)){
			newick = line;
		}
	}
	//Reads in sequences ok
	sequence_length = 1000;
	tree = TreeUtil.readTree(newick);
	SGen.TakeInTree(tree, sequence_length);

}
