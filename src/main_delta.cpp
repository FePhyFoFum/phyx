/*
 * main_delta.cpp
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "delta.h"
#include "log.h"

#include "omp.h"

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    TreeReader tr;

    if (argc != 4 && argc != 3) {
        cout << "usage: phyx_delta l r o" << endl;
        cout << "    OR" << endl;
        cout << "usage: phyx_delta newickfile outfile" << endl;
        exit(0);
    }
    if (argc == 4) {
        Delta delta;
        vector<double> nums = delta.delta(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
        cout << nums[0] << " " << nums[1] << endl;
    }
    if (argc == 3) {
        TreeReader tr;
        Tree * tree;
        ifstream ifs( argv[1] );
        string temp;
        int count = 1;
        while ( getline( ifs, temp ) ) {
            if (temp.size() > 1) {
                tree = tr.readTree(temp);
                cout << "Tree "<< count <<" has " << tree->getExternalNodeCount()
                    << " leaves." << endl;
                //ret.push_back(intree);
                count++;
            }
        }
        //moving through the tree

        Delta delta;
        #pragma omp parallel for ordered num_threads(4)
        for (int i=0; i < tree->getInternalNodeCount(); i++) {
            if (tree->getInternalNode(i)->get_num_leaves() > 1
                    && tree->getInternalNode(i)->isRoot() == false
                    && tree->getInternalNode(i)->getChildCount() == 2) {
                int left = tree->getInternalNode(i)->getChild(0)->get_num_leaves();
                int right = tree->getInternalNode(i)->getChild(1)->get_num_leaves();
                int out = tree->getInternalNode(i)->getParent()->get_num_leaves()-(left+right);
                vector <double> nums = delta.delta(left, right, out);
                cout << "left:" << left << " right:" << right << " out:"
                    << out << " delta:" << nums[0] << " p:" << nums[1] << endl;
                std::stringstream sout;
                sout << nums[1];
                tree->getInternalNode(i)->setName(sout.str());
            }
        }
        ofstream outTreeFile;
        outTreeFile.open(argv[2],ios::app );
        outTreeFile << tree->getRoot()->getNewick(true) << ";" << endl;
        outTreeFile.close();
    }
    return EXIT_SUCCESS;
}
