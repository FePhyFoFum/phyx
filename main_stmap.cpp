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
#include "optimize_state_reconstructor.h"
#include "optimize_nlopt.h"

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

	int nstates  = seqs[0].get_sequence().length();
	cout << "states: " << nstates << endl;
	RateModel rm(nstates);
	rm.setup_P(0.1,false);
	StateReconstructor sr(rm);
	sr.set_tree(tree);
	sr.set_tip_conditionals(seqs);

	mat free_var(nstates,nstates);free_var.fill(0);
	int ct = 0;
	for(int i=0;i<nstates;i++){
		for(int j=0;j<nstates;j++){
			if(i!=j){
				free_var(i,j) = ct;
				ct += 1;
			}
		}
	}
	cout << free_var << endl;
	cout << ct << endl;

//	OptimizeStateReconstructor osr(&rm,&sr,&free_var,1);
//	free_var = osr.optimize();
//	rm.setup_Q(free_var);
//	cout << free_var << endl;

	optimize_sr_nlopt(&rm,&sr,&free_var,ct);
	cout << free_var << endl;
	rm.setup_Q(free_var);

	sr.set_store_p_matrices(true);
	cout << sr.eval_likelihood() << endl;

/*	sr.prepare_ancstate_reverse();
	vector<double> lhoods;
	for(int i=0;i<tree->getInternalNodeCount();i++){
		cout <<"node: " << tree->getInternalNode(i)->getName() << endl;
		lhoods = sr.calculate_ancstate_reverse(*tree->getInternalNode(i));
		cout << lhoods[0]/calculate_vector_double_sum(lhoods) << " ";
		cout << lhoods[1]/calculate_vector_double_sum(lhoods) << " ";
		cout << lhoods[2]/calculate_vector_double_sum(lhoods) << " ";
		cout << lhoods[3]/calculate_vector_double_sum(lhoods) << endl;
	}


	cout << tree->getInternalNode(3)->getName() << endl;
	sr.prepare_stochmap_reverse_all_nodes(0,0);
	sr.prepare_ancstate_reverse();
	vector<double> stoch = sr.calculate_reverse_stochmap(*tree->getInternalNode(3),true);
	cout << calculate_vector_double_sum(stoch)/calculate_vector_double_sum(lhoods) << endl;
	sr.prepare_stochmap_reverse_all_nodes(1,2);
	sr.prepare_ancstate_reverse();
	stoch = sr.calculate_reverse_stochmap(*tree->getInternalNode(3),true);
	cout << calculate_vector_double_sum(stoch)/calculate_vector_double_sum(lhoods) << endl;
	sr.prepare_stochmap_reverse_all_nodes(2,3);
	sr.prepare_ancstate_reverse();
	stoch = sr.calculate_reverse_stochmap(*tree->getInternalNode(3),true);
	cout << calculate_vector_double_sum(stoch)/calculate_vector_double_sum(lhoods) << endl;
	sr.prepare_stochmap_reverse_all_nodes(3,0);
	sr.prepare_ancstate_reverse();
	stoch = sr.calculate_reverse_stochmap(*tree->getInternalNode(3),true);
	cout << calculate_vector_double_sum(stoch)/calculate_vector_double_sum(lhoods) << endl;
*/
	return EXIT_SUCCESS;
}
