/*
 * main_upgma.cpp
 *
 *  Created on: Jun 10, 2015
 *      Author: joe
 */



#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <map>
#include <iterator>
#include "upgma.h"
#include "node.h"
#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "concat.h"
using namespace std;

void print_help() {
    cout << "Basic UPGMA Tree Maker." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << "Individual files can be of different formats." << endl;
    cout << endl;
    cout << "Usage: pxupgma [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     list of input sequence files" << endl;
    cout << " -o, --outf=FILE     output newick file, stout otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxconcat 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

/*
static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};
*/
//int main(int argc, char * argv[]) {
int main(){

	map<string, string> sequences;
	map<string, string>::iterator iter;
	map<int, string> NameKey;
	vector< vector<double> > Matrix;
	vector<string> names;
	string fasta;
	int count = 0;
	//fasta = ("TestFiles/drosophila.aln");
	fasta = ("TestFiles/Real_Test.fa");
	//cin >> fasta;
	UPGMA functions;
	sequences = functions.FastaToOneLine(fasta);
	for(iter = sequences.begin(); iter != sequences.end(); iter++){
		NameKey[count] = iter -> first;
		names.push_back(iter -> first);
		count++;
	}
	Matrix = functions.BuildMatrix(sequences);
	functions.TREEMAKE(names, NameKey, Matrix);
	return 0;
}
