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
#include <getopt.h>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "bd_sim.h"

/*
void print_help(){
    cout << "Get information about an mrca" << endl;
    cout << endl;
    cout << "Usage: pxmrca [OPTION]... [FILE]..."<<endl;
    cout << endl; 
    cout << " -m, --mean=VALUE    mean value under which seqs are filtered" << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" <<endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>"<<endl;
}

string versionline("pxmrca 0.1\nCopyright (C) 2013 FePhyFoFum\nLiscence GPLv2\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"mean", required_argument, NULL, 'm'},
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

*/

int main(int argc, char * argv[]){
    TreeReader tr;
    
    if (argc != 3){
        cout << "usage: pxmrca MRCA newickfile" << endl;
        exit(0);
    }
    
    ifstream infile(argv[1]);
    if (!infile){
        cerr << "Could not open mrcafile." << endl;
        return 1;
    }
    
    string mrcaline;
    map<string, vector<string> > mrcas;
    while (getline(infile, mrcaline)){
        vector<string> searchtokens;
        tokenize(mrcaline, searchtokens, "     ");
        cout << "Read in " << searchtokens.size() << " tokens!" << endl;
        for(unsigned int j=0; j < searchtokens.size(); j++){
            trim_spaces(searchtokens[j]);
        }
        vector<string> vec;
        vec.push_back(searchtokens[1]);
        vec.push_back(searchtokens[2]);
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
        cout << (*it).first<<" "<< nd->get_num_leaves() << " " << nd->getName() << endl;
    }

    delete tree;
    return EXIT_SUCCESS;
}
