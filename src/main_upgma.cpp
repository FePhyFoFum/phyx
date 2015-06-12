/*
 * main_upgma.cpp
 *
 *  Created on: Jun 10, 2015
 *      Author: joe
 */


// TODO: need to remove unnecessary includes
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
//#include <algorithm>
#include <map>
#include <iterator>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "upgma.h"
//#include "node.h"
//#include "sequence.h"
//#include "seq_reader.h"
#include "utils.h"
//#include "concat.h"

void print_help() {
    cout << "Basic UPGMA Tree Maker." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxupgma [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output newick file, stout otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxupgma 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");


static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    bool fileset = false;
    bool outfileset = false;
    string seqf = "";
    string outf = ""; // not used at the moment
    
    while (1) {
        int oi = -1;
        int curind = optind;
        int c = getopt_long(argc, argv, "s:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
            case 'V':
                cout << versionline << endl;
                exit(0);
            default:
                print_error(argv[0],(char)c);
                exit(0);
        }
    }
    
    // only taking files at the moment (not stdin)
    if (!fileset) {
        cout << "you must specify an input sequence file" << endl;
        exit(0);
    }
    
    map<string, string> sequences;
    map<string, string>::iterator iter;
    map<int, string> NameKey;
    vector< vector<double> > Matrix;
    vector<string> names;
    string fasta = seqf; // temporary
    int count = 0;
    //fasta = ("TestFiles/drosophila.aln");
    //fasta = ("TestFiles/Real_Test.fa");
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
    cout << "Newick:" << endl << functions.get_newick() << endl;
    
    cout << endl;
    
    // alternate:
    UPGMA terp(seqf);
    cout << "Newick:" << endl << terp.get_newick() << endl;
    
    //return 0;
    return EXIT_SUCCESS;
}
