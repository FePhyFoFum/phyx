/*
 * main_nni.cpp
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
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

void print_help() {
    cout << "Nearest Neighbor Interchange Program" << endl;
    cout << "This will take newick of nexus files" << endl;
    cout << endl;
    cout << "Usage: pxnni [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE    input tree file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output tree file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]){

    bool outfileset = false;
    bool fileset = false;
    char * outf;
    char * seqf;
	/*
    if (argc > 2){
        cout << "usage: pxnni newickfile" << endl;
        exit(0);
    }*/
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
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
                cout << "Version" << endl;
                exit(0);
            default:
                //print_error(argv[0],(char)c);
                cout << "? Maybe try the help menu -h" << endl;
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
        fstr = new ifstream(seqf);
        pios = fstr;
    } else {
        pios = &cin;
    }
    srand(time(0));
    TreeReader tr;
    vector<string> lines;

    //reading from standard input for piping
    /*
    if(argc == 1){
        for (std::string line; std::getline(std::cin, line);) {
            lines.push_back(line);
        }
    }*/
    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "this really only works with nexus or newick" << endl;
        exit(0);
    }
    //reading from a file
    /*
    if (argc == 2){
        ifstream infile(argv[1]);
        if (!infile){
            cerr << "Could not open treefile." << endl;
            return 1;
        }
        string line;
        while (getline(infile, line)){
            lines.push_back(line);
        }
        infile.close();
    }*/
    /*while (getline(infile, fstr)){
        lines.push_back(line);
    }
    infile.close();*/
    int treeCounter = 0;
    bool going = true;
    if (ft == 1) { // newick. easy
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
            if (tree != NULL) {
                //cout << "Working on tree #" << treeCounter << endl;
                map<Node*,vector<Node*> > tree_map;
				create_tree_map_from_rootnode(tree,tree_map);
				nni_from_tree_map(tree,tree_map);
				(*poos) << tree->getRoot()->getNewick(true) << endl;
                delete tree;
                treeCounter++;
            }
        }
    }
    /*
    Tree * tree = tr.readTree(lines[0]);
    map<Node*,vector<Node*> > tree_map;
    create_tree_map_from_rootnode(tree,tree_map);
    nni_from_tree_map(tree,tree_map);
	*/
    //cout << tree->getRoot()->getNewick(true) << ";" << endl;
    //delete tree;
    return EXIT_SUCCESS;
}
