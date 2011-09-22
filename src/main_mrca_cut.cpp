/*
 * main_mrca_cut.cpp
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

	if (argc != 3){
		cout << "usage: phyx_mrca_cut MRCA newickfile" << endl;
		exit(0);
	}

	ifstream infile(argv[1]);
	if (!infile){
		cerr << "Could not open mrcafile." << endl;
		return 1;
	}

	string mrcaline;
	map<string,vector<string> > mrcas;
	while (getline(infile, mrcaline)){
		vector<string> searchtokens;
		Tokenize(mrcaline, searchtokens, " 	");
		for(unsigned int j=0;j<searchtokens.size();j++){
			TrimSpaces(searchtokens[j]);
		}
		vector<string> vec;vec.push_back(searchtokens[1]);vec.push_back(searchtokens[2]);
		mrcas[searchtokens[0]] = vec;
	}
	infile.close();



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
	cout << tree->getExternalNodeCount() << endl;

	map<string,vector<string> >::iterator it;
	for (it=mrcas.begin(); it!=mrcas.end(); it++){
		Node * nd = tree->getMRCA((*it).second);
		ofstream outFile;
		outFile.open((*it).first.c_str(), ios::out);
		outFile << nd->getNewick(true) << endl;
		outFile.close();
	}

	delete tree;
	return EXIT_SUCCESS;
}
