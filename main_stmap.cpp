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
#include "sequence.h"
#include "phylip_reader.h"
#include "state_reconstructor.h"
#include "rate_model.h"

int main(int argc, char * argv[]){
	TreeReader tr;

	if (argc != 3){
		cout << "usage: phyx_stmap seqfile newickfile" << endl;
		exit(0);
	}

	vector<Sequence> seqs;
	PhylipReader pr;
	bool phyl = pr.readFile(argv[1],seqs);
	if(phyl == false){
		cout << "the sequence file is not phylip" << endl;
		exit(0);
	}
	cout << "sequences: " << seqs.size() << endl;



	ifstream infile2(argv[2]);
	if (!infile2){
		cerr << "Could not open treefile." << endl;
		return 1;
	}
	vector<string> lines;
	string line;
	while (getline(infile2, line)){
		lines.push_back(line);
	}
	infile2.close();

	Tree * tree = tr.readTree(lines[0]);
	cout << "tips: "<< tree->getExternalNodeCount() << endl;

	cout << "states: " << seqs[0].get_sequence().length() << endl;
	RateModel rm(seqs[0].get_sequence().length());
	rm.setup_P(0.1,false);
	StateReconstructor sr(rm);


	return EXIT_SUCCESS;
}
