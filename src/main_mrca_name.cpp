/*
 * main_mrca_name.cpp
 *
 */

/*
 * The idea behind this is to allow for the naming of internal nodes based
 * on given MRCAS and a set of names in a file the input of which should
 * look like
 * mrca list seperated by spaces\tname
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

    if (argc != 4){
	cout << "usage: phyx_mrca_names MRCA newickfile outfile" << endl;
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
	tokenize(mrcaline, searchtokens, "\t");
	for(unsigned int j=0;j<searchtokens.size();j++){
	    trim_spaces(searchtokens[j]);
	}
	//first searchtoken is the set of names
	vector<string> vec;
	vector<string> searchtokens2;
	tokenize(searchtokens[0],searchtokens2," ");
	for(unsigned int j=0;j<searchtokens2.size();j++){
	    trim_spaces(searchtokens2[j]);
	    vec.push_back(searchtokens2[j]);
	}
//second searchtoken is the set of mrca names
	mrcas[searchtokens[1]] = vec;
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
	nd->setName((*it).first);
    }
    ofstream outtreefile;
    outtreefile.open(argv[3],ios::app);
    outtreefile << tree->getRoot()->getNewick(true) << ";" << endl;
    outtreefile.close();
    delete tree;
    return EXIT_SUCCESS;
}
