#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <nlopt.hpp>

using namespace std;

#include "ls_tr.h"
#include "tree.h"
#include "node.h"
#include "tree_utils.h"


TreeInfo::TreeInfo (Tree * intree) {
    tree = intree;
    calc_stats();
}

TreeInfo::TreeInfo (Tree * intree, bool const& ultracheck, bool const& binarycheck,
    bool const& agecheck, bool const& rootedcheck, bool const& ntipcheck,
    bool const& lengthcheck, ostream* poos) {
    tree = intree;
    if (ultracheck) {
        ultrametric_tree = is_ultrametric_paths(tree);
        (*poos) << std::boolalpha << ultrametric_tree << endl;
    } else if (binarycheck) {
        binary_tree = is_binary(tree);
        (*poos) << std::boolalpha << binary_tree << endl;
    } else if (agecheck) {
        ultrametric_tree = is_ultrametric_paths(tree);
        if (ultrametric_tree) {
            rootheight = tree->getRoot()->getHeight();
            (*poos) << rootheight << endl;
        } else {
            (*poos) << "NA" << endl;
        }
    } else if (rootedcheck) {
        rooted_tree = is_rooted(tree);
        (*poos) << std::boolalpha << rooted_tree << endl;
    } else if (ntipcheck) {
        ntips = tree->getExternalNodeCount();
        (*poos) << ntips << endl;
    } else if (lengthcheck) {
        treelength = get_tree_length(tree);
        has_branchlengths = (treelength > 0.0) ? true : false;
        if (has_branchlengths) {
            (*poos) << treelength << endl;
        } else {
            (*poos) << "NA" << endl;
        }
    } 
}

void TreeInfo::calc_stats () {
    treelength = get_tree_length(tree);
    has_branchlengths = (treelength > 0.0) ? true : false;
    nintnodes = tree->getInternalNodeCount();
    ntips = tree->getExternalNodeCount();
    rooted_tree = is_rooted(tree);
    binary_tree = is_binary(tree);
    ultrametric_tree = is_ultrametric_paths(tree);
    if (ultrametric_tree) {
        rootheight = tree->getRoot()->getHeight();
    } else {
        rootheight = 0.0;
    }
}


void TreeInfo::get_stats (ostream* poos) {
    (*poos) << "rooted: " << std::boolalpha << rooted_tree << endl;
    (*poos) << "binary: " << std::boolalpha << binary_tree << endl;
    (*poos) << "nterminal: " << ntips << endl;
    (*poos) << "ninternal: " << nintnodes << endl;
    (*poos) << "branch lengths: " << std::boolalpha << has_branchlengths << endl;
    if (has_branchlengths) {
        (*poos) << "treelength: " << treelength << endl;
        (*poos) << "ultrametric: " << std::boolalpha << ultrametric_tree << endl;
    } else {
        (*poos) << "treelength: NA" << endl;
        (*poos) << "ultrametric: NA" << endl;
    }
    
    if (ultrametric_tree) {
        (*poos) << "rootheight: " << rootheight << endl;
    } else {
        (*poos) << "rootheight: NA" << endl;
    }
}
