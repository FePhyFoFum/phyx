/*
 * main_nni.cpp
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "tree_utils.h"

int main(int argc, char * argv[]){

    if (argc > 2){
        cout << "usage: pxnni newickfile" << endl;
        exit(0);
    }

    srand(time(0));
    TreeReader tr;
    vector<string> lines;

    //reading from standard input for piping
    if(argc == 1){
        for (std::string line; std::getline(std::cin, line);) {
            lines.push_back(line);
        }
    }
    //reading from a file
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
    }

    Tree * tree = tr.readTree(lines[0]);
    map<Node*,vector<Node*> > tree_map;
    create_tree_map_from_rootnode(tree,tree_map);
    nni_from_tree_map(tree,tree_map);

    cout << tree->getRoot()->getNewick(true) << ";" << endl;
    delete tree;
    return EXIT_SUCCESS;
}
