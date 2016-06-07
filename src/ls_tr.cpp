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
    (*poos) << "treelength: " << treelength << endl;
    (*poos) << "ultrametric: " << std::boolalpha << ultrametric_tree << endl;
    if (ultrametric_tree) {
        (*poos) << "rootheight: " << rootheight << endl;
    } else {
        (*poos) << "rootheight: NA" << endl;
    }
}
