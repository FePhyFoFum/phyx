// g++ -std=c++11 main_JWB.cpp tree_reader.cpp utils.cpp node.cpp superdouble.cpp tree.cpp -o JWB

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "tree_reader.h"
#include "tree.h"
#include "utils.h"
#include "node.h"
#include "superdouble.h"

int main(int argc, char * argv[]) {

    if (argc != 2) {
        cout << "usage: ./JWB treefile" << endl;
        exit(0);
    }
    
    char * treef = argv[1];
    
    istream * pios;
    ostream * poos;
    ifstream * fstr;
    ofstream * ofstr;
    
    poos = &cout;
    fstr = new ifstream(treef);
    pios = fstr;
    
    TreeReader tr;
    string retstring;
    
    int ft = test_tree_filetype_stream(*pios, retstring);
    
    int treeCounter = 0;
    bool going = true;
    if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
            if (tree != NULL) {
                (*poos) << "original tree: " << endl;
                (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
                
                Tree * atree = new Tree(*tree);
                Tree foo;
                *foo = *tree;
                atree = tree;
                //atree->unRoot();
                foo.unRoot();
                (*poos) << "altered tree tree: " << endl;
                (*poos) << atree->getRoot()->getNewick(true) << ";" << endl;
                
                (*poos) << "and the original tree again: " << endl;
                (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
                
                (*poos) << "foo: " << endl;
                (*poos) << foo.getRoot()->getNewick(true) << ";" << endl;
                //delete tree;
                //delete atree;
                treeCounter++;
            }
        }
    }
    
    
    
    return EXIT_SUCCESS;
}
