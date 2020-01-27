// this is not present in the Makefile


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

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
        std::cout << "usage: phyx_delta l r o" << std::endl;
        std::cout << "    OR" << std::endl;
        std::cout << "usage: phyx_delta newickfile outfile" << std::endl;
        exit(0);
    }
    if (argc == 4) {
        Delta delta;
        vector<double> nums = delta.delta(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
        std::cout << nums[0] << " " << nums[1] << std::endl;
    }
    if (argc == 3) {
        TreeReader tr;
        Tree * tree;
        ifstream ifs( argv[1] );
        std::string temp;
        int count = 1;
        while ( getline( ifs, temp ) ) {
            if (temp.size() > 1) {
                tree = tr.readTree(temp);
                std::cout << "Tree "<< count <<" has " << tree->getExternalNodeCount()
                    << " leaves." << std::endl;
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
                std::vector<double> nums = delta.delta(left, right, out);
                std::cout << "left:" << left << " right:" << right << " out:"
                    << out << " delta:" << nums[0] << " p:" << nums[1] << std::endl;
                std::stringstream sout;
                sout << nums[1];
                tree->getInternalNode(i)->setName(sout.str());
            }
        }
        ofstream outTreeFile;
        outTreeFile.open(argv[2],ios::app );
        outTreeFile << tree->getRoot()->getNewick(true) << ";" << std::endl;
        outTreeFile.close();
    }
    return EXIT_SUCCESS;
}
