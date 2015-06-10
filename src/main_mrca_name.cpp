/*
 * main_mrca_name.cpp
 *
*/

/*
 * The idea behind this is to allow for the naming of internal nodes based
 * on given MRCAS and a set of names in a file the input of which should
 * look like:
 * MRCANAME = tip1 tip2 ...
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "tree_utils.h"

void print_help(){
    cout << "Label internal nodes with clade names." << endl;
    cout << "Takes in newick tree and MRCA file with format:" << endl;
    cout << "MRCANAME = tip1 tip2 ..." << endl;
    cout << endl;
    cout << "Usage: pxmrcaname [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE    input newick tree file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output newick file, stout otherwise" << endl;
    cout << " -m, --mrca=FILE     file containing MRCA declarations" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxmrcaname 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv2\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"mrca", required_argument, NULL, 'm'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};


int main(int argc, char * argv[]){
    TreeReader tr;
    char * outf;
    char * treef;
    char * mrcaf;
    bool outfileset = false;
    bool fileset = false;
    bool mrcaset = false;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:m:hV", long_options, &oi);
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
            case 'm':
                mrcaset = true;
                mrcaf = strdup(optarg);
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
    
    if (!mrcaset) {
        cout << "Must supply mrca file." << endl;
        exit(0);
    }
    
    istream* pios;
    ostream* poos;
    ifstream* fstr;
    ofstream* ofstr;
    
    if (outfileset == true) {
        ofstr = new ofstream(outf, ios::app);
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
    
    /* 
       collect clade names
       expecting (new) format:
       MRCANAME = tip1 tip2 ... 
    */
    ifstream inmrca(mrcaf);
    string mrcaline;
    map<string, vector<string> > mrcas;
    while (getline(inmrca, mrcaline)){
        if (mrcaline.empty()) {
            continue;
        }
        vector<string> searchtokens;
        tokenize(mrcaline, searchtokens, "=");
        string mrcaname = searchtokens[0];
        trim_spaces(mrcaname);
        searchtokens.erase(searchtokens.begin());
        searchtokens = tokenize(searchtokens[0]);
        mrcas[mrcaname] = searchtokens;
    }
    inmrca.close();
    
// collect tree(s)
    vector<string> lines;
    string line;
    while (getline(*pios, line)){
        lines.push_back(line);
    }
    
// allow multiple trees
    for (unsigned int i = 0; i < lines.size(); i++) {
        Tree * tree = tr.readTree(lines[i]);
        map<string, vector<string> >::iterator it;
        for (it = mrcas.begin(); it != mrcas.end(); it++){
            //cout << "Dealing with clade '" << (*it).first << "'" << endl;
            if (!check_names_against_tree(tree, (*it).second)) {
                cout << "Check mrca file for typos." << endl;
                exit (0);
            }
            Node * nd = tree->getMRCA((*it).second);
            nd->setName((*it).first);
        }
        (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
        delete tree;
    }
    
    if (fileset) {
        fstr->close();
        delete pios;
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
