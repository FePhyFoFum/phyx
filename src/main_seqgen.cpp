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

//Standalone compile line, needs armadillo 5.2
//g++ -std=c++11 branch_segment.cpp seqgen.cpp node.cpp tree.cpp tree_reader.cpp main_seqgen.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp -larmadillo  -llapack -lblas -lpthread -lm -o test
void print_help() {
    cout << "Basic Sequence Simulator under the GTR Model." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxseqsim [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE     input treefile, stdin otherwise" << endl;
    cout << " -c, --comp=Input     input base composition comma separated in order A,T,C, Default is .25 for all" << endl;
    cout << " -a, --ancestor=Input  Prints the ancestral sequences at each node, Default if off" << endl;
    cout << " -l, --len=Input     length of the sequences that you would like to generate" << endl;
    cout << " -r, --ratemat=Input  input values for rate matrix to evolve under, Default is JC69" << endl;
    cout << " -o, --outf=FILE     output newick file, stout otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxseqsim 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");


static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
	{"comp", required_argument, NULL, 'c'},
	{"length", required_argument, NULL, 'l'},
	{"Ancestors", required_argument, NULL, 'l'},
	{"ratemat", required_argument, NULL, 'r'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(){


	TreeReader TreeUtil;
	Tree * tree;
	SEQGEN SGen;
	double length;
	string path, newick, line;
	int sequence_length = 0;
	/*
	 * Default Base Frequencies and Rate Matrix
	 *
	 */
	vector<double> basefreq(4, 0.0);
    basefreq[0] = .25;
    basefreq[1] = .25;
    basefreq[2] = .25;
    basefreq[3] = (((1 - basefreq[0]) - basefreq[1]) - basefreq[2]);
    vector< vector<double> > rmatrix(4, vector<double>(4, 1.0));
    for (int i = 0; i < rmatrix.size(); i++){
    	for (int j = 0; j < rmatrix.size(); j++){

    		if (i == j){//Fill Diagnol

    			rmatrix[i][j] *= -1.0;

    		}else{
    			rmatrix[i][j] *= 0.33333;
    		}
    	}
    }
	ifstream readline;
	path = ("TestFiles/RAxML_bestTree.drosML");
	readline.open(path.c_str());
	if (readline.is_open()){
		while (getline(readline, line)){
			newick = line;
		}
	}
	sequence_length = 1001;
	tree = TreeUtil.readTree(newick);
	SGen.TakeInTree(rmatrix, tree, sequence_length, basefreq);

}
