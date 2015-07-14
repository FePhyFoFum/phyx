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
#include <cstring>
#include <getopt.h>

using namespace std;

#include "seqgen.h"
//#include "sequence.h"
//#include "seq_reader.h"
#include "utils.h"
//#include "node.h"
#include "tree.h"
#include "tree_reader.h"

//Standalone compile line, needs armadillo 5.2
//g++ -std=c++11 branch_segment.cpp seqgen.cpp node.cpp tree.cpp tree_reader.cpp main_seqgen.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp -larmadillo  -llapack -lblas -lpthread -lm -o test

// example call: ./pxseqgen -t TEST/ultra_100.tre -l 25


void print_help() {
    cout << "Basic Sequence Simulator under the GTR Model." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxseqgen [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE     input treefile, stdin otherwise" << endl;
    cout << " -c, --comp=Input     comma-delimited base freqs in order: A,T,C,G. default is equal" << endl;
    cout << " -a, --ancestors      print the ancestral sequences at each node. default is no" << endl;
    cout << " -l, --len=INT        length of sequences to generate. default is 1000" << endl;
    cout << " -r, --ratemat=Input  input values for rate matrix. default is JC69" << endl;
    cout << " -o, --outf=FILE      output seq file, stout otherwise" << endl;
    cout << " -n, --nreps=INT      number of replicates" << endl;
    cout << " -x, --seed=INT       random number seed, clock otherwise" << endl;
    cout << "     --help           display this help and exit" << endl;
    cout << "     --version        display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxseqgen 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"comp", required_argument, NULL, 'c'},
    {"length", required_argument, NULL, 'l'},
    {"ancestors", no_argument, NULL, 'a'},
    {"ratemat", required_argument, NULL, 'r'},
    {"nreps", required_argument, NULL, 'n'},
    {"seed", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    bool outfileset = false;
    bool fileset = false;
    bool showancs = false;
    int seqlen = 1000;
    string infreqs;
    string inrates;
    char * outf;
    char * treef;
    vector <double> basefreq(4, 0.25);
    vector <double> userrates(4, 0.25);
    int nreps = 1; // not implemented at the moment
    int seed = -1;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:c:l:ar:n:x:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'c':
                infreqs = strdup(optarg);
                parse_double_list(infreqs, basefreq, ", ");
                break;
            case 'l':
                seqlen = atoi(strdup(optarg));
                break;
            case 'a':
                showancs = true;
                break;
            case 'r':
                inrates = strdup(optarg);
                parse_double_list(inrates, userrates, ", ");
                break;
            case 'n':
                nreps = atoi(strdup(optarg));
                break;
            case 'x':
                seed = atoi(strdup(optarg));
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
    
    istream* pios;
    ostream* poos;
    ifstream* fstr;
    ofstream* ofstr;
    
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    if (fileset == true) {
        fstr = new ifstream(treef);
        pios = fstr;
    } else {
        pios = &cin;
    }
    
    //TreeReader TreeUtil;
    //Tree * tree;
    //SEQGEN SGen(seqlen, basefreq);
    //double length;
    string path, newick, line;
    //int sequence_length = 0;
    /*
     * Default Base Frequencies and Rate Matrix
     *
     */
    //vector <double> basefreq(4, 0.0);
    //basefreq[0] = .25;
    //basefreq[1] = .25;
    //basefreq[2] = .25;
    //basefreq[3] = 1.0 - basefreq[0] - basefreq[1] - basefreq[2];
    vector< vector <double> > rmatrix(4, vector<double>(4, 0.33));
    for (unsigned int i = 0; i < rmatrix.size(); i++) {
        for (unsigned int j = 0; j < rmatrix.size(); j++) {
            if (i == j) {//Fill Diagnol
                rmatrix[i][j] = -1.0;
            }
        }
    }
//    ifstream readline;
//    //path = ("TestFiles/RAxML_bestTree.drosML");
//    path = ("TEST/ultra_100.tre");
//    readline.open(path.c_str());
//    if (readline.is_open()) {
//        while (getline(readline, line)) {
//            newick = line;
//        }
//    }
//    //sequence_length = 1001;
//    tree = TreeUtil.readTree(newick);
    
    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "this really only works with nexus or newick" << endl;
        exit(0);
    }
    
    // allow > 1 tree in input. passing but not yet using nreps
    int treeCounter = 0;
    bool going = true;
    bool exists;
    if (ft == 1) { // newick. easy
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
            if (tree != NULL) {
                cout << "Working on tree #" << treeCounter << endl;
                SEQGEN SGen(seqlen, basefreq, rmatrix, tree, nreps);
                delete tree;
                treeCounter++;
            }
        }
    } else if (ft == 0) { // Nexus. need to worry about possible translation tables
        map <string, string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (going == true) {
                cout << "Working on tree #" << treeCounter << endl;
                SEQGEN SGen(seqlen, basefreq, rmatrix, tree, nreps);
                delete tree;
                treeCounter++;
            }
        }
    }
    
    
    
//    while (going) {
//        Tree * tree;
//        if (ft == 1) { // newick
//            tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
//            if (tree != NULL) {
//                cout << "Working on tree #" << treeCounter << endl;
//                SEQGEN SGen(seqlen, basefreq, rmatrix, tree, nreps);
//                delete tree;
//                treeCounter++;
//            }
//        } else { // Nexus
//            vector<string> retstrings;
//            retstrings.push_back(retstring);
//            map<string,string> translation_table;
//            bool ttexists;
//            ttexists = get_nexus_translation_table(*pios, &translation_table, &retstrings);
//            if (retstrings.size() > 0) {
//                retstring = retstrings[retstrings.size()-1];
//            }
//        }
//    }
    
    //SEQGEN SGen(seqlen, basefreq, rmatrix);
    //SGen.TakeInTree(rmatrix, tree);
    return EXIT_SUCCESS;
}
