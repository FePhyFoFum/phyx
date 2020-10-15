#include "tree.h"
#include "node.h"
#include "tree_utils.h"
#include "clean_tree.h"


CleanTree::CleanTree (Tree * intree, bool remove_root_edge, bool remove_labels,
        bool remove_knuckles) {
    tree_ = intree;
    root_ = tree_->getRoot();
    remove_root_edge_ = remove_root_edge;
    remove_labels_ = remove_labels;
    remove_knuckles_ = remove_knuckles;
    clean_properties();
}

// remove annotations, 'knuckles', and root edges
void CleanTree::clean_properties () {
    if (remove_labels_) {
        if (tree_->hasNodeAnnotations()) {
            //cout << "Found node annotations!" << endl;
            remove_annotations(tree_);
        }
        if (tree_->hasNodeNames()) {
            //cout << "Found node names!" << endl;
            remove_internal_names(tree_);
        }
    }
    if (remove_knuckles_) {
        deknuckle_tree(tree_);
    }
    if (remove_root_edge_) {
        if (has_root_edge(tree_)) {
            //cout << "Tree has a root edge." << endl;
            tree_->removeRootEdge();
        }
    }
}
