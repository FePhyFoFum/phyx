/*
 * main_strec.cpp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "sequence.h"
#include "fasta_util.h"

int main(int argc, char * argv[]){
	TreeReader tr;

	if (argc != 1){
		cout << "usage: phyx_test" << endl;
		exit(0);
	}
	string datafile = "../../../projects/PHLAWD_fish/12S.keep";
	vector<Sequence> seqs;
	FastaUtil pr;
	bool phyl = pr.readFile(datafile,seqs);
	cout << "sequences: " << seqs.size() << endl;
	cout << seqs[0].get_sequence() <<endl;
	seqs[0].reverse_complement();
	cout << endl;
	cout << seqs[0].get_sequence()<< endl;
	cout << "writing file" << endl;
	string outfile = "test.fasta";
	pr.writeFileFromVector(outfile,seqs);
	return EXIT_SUCCESS;
}
