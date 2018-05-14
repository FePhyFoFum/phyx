#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "tree.h"
#include "node.h"
#include "tree_utils.h"
#include "clean_tree.h"

CleanTree::CleanTree (Tree * intree) {
    tree_ = intree;
    root_ = tree_->getRoot();
    clean_all();
}

// remove annotations, 'knuckles', and root edges
void CleanTree::clean_all () {
    if (tree_->hasNodeAnnotations()) {
        //cout << "Found node annotations!" << endl;
        remove_annotations(tree_);
    }
    if (tree_->hasNodeNames()) {
        //cout << "Found node names!" << endl;
        remove_internal_names(tree_);
    }
    deknuckle_tree(tree_);
    if (has_root_edge(tree_)) {
        //cout << "Tree has a root edge." << endl;
        tree_->removeRootEdge();
    }
}