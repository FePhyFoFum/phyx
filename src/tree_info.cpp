#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "tree_info.h"
#include "tree.h"
#include "node.h"
#include "tree_utils.h"


TreeInfo::TreeInfo (Tree * intree) {
    tree_ = intree;
    calc_stats();
}


TreeInfo::TreeInfo (Tree * intree, bool const& ultracheck, bool const& binarycheck,
        bool const& agecheck, bool const& rootedcheck, bool const& ntipcheck,
        bool const& lengthcheck, bool const& namecheck, bool const& rtvarcheck,
        ostream* poos) {
    tree_ = intree;
    if (ultracheck) {
        ultrametric_tree_ = is_ultrametric_paths(tree_);
        (*poos) << std::boolalpha << ultrametric_tree_ << endl;
    } else if (rtvarcheck) {
        has_branchlengths_ = tree_->hasEdgeLengths();
        if (has_branchlengths_) {
            rtvar_ = get_root_tip_var(tree_);
            (*poos) << rtvar_ << endl;
        } else {
            (*poos) << "NA" << endl;
        }
    } else if (binarycheck) {
        binary_tree_ = is_binary(tree_);
        (*poos) << std::boolalpha << binary_tree_ << endl;
    } else if (agecheck) {
        ultrametric_tree_ = is_ultrametric_paths(tree_);
        if (ultrametric_tree_) {
            rootheight_ = tree_->getRoot()->getHeight();
            (*poos) << rootheight_ << endl;
        } else {
            (*poos) << "NA" << endl;
        }
    } else if (rootedcheck) {
        rooted_tree_ = is_rooted(tree_);
        (*poos) << std::boolalpha << rooted_tree_ << endl;
    } else if (ntipcheck) {
        ntips_ = tree_->getExternalNodeCount();
        (*poos) << ntips_ << endl;
    } else if (lengthcheck) {
        has_branchlengths_ = tree_->hasEdgeLengths();
        if (has_branchlengths_) {
            treelength_ = get_tree_length(tree_);
            (*poos) << treelength_ << endl;
        } else {
            (*poos) << "NA" << endl;
        }
    } else if (namecheck) {
        tip_labels_ = get_tip_labels(tree_);
        for (unsigned int i = 0; i < tip_labels_.size(); i++) {
            (*poos) << tip_labels_[i] << endl;
        }
    }
}


void TreeInfo::calc_stats () {
    has_branchlengths_ = tree_->hasEdgeLengths();
    if (has_branchlengths_) {
        treelength_ = get_tree_length(tree_);
    } else {
        treelength_ = 0.0;
    }
    nintnodes_ = tree_->getInternalNodeCount();
    ntips_ = tree_->getExternalNodeCount();
    rooted_tree_ = is_rooted(tree_);
    binary_tree_ = is_binary(tree_);
    ultrametric_tree_ = is_ultrametric_paths(tree_);
    if (rooted_tree_){
        rtvar_ = get_root_tip_var(tree_);
    } else {
        rtvar_ = 0.0;
    }
    if (ultrametric_tree_) {
        rootheight_ = tree_->getRoot()->getHeight();
    } else {
        rootheight_ = 0.0;
    }
}


void TreeInfo::get_stats (ostream* poos) {
    (*poos) << "rooted: " << std::boolalpha << rooted_tree_ << endl;
    (*poos) << "binary: " << std::boolalpha << binary_tree_ << endl;
    (*poos) << "nterminal: " << ntips_ << endl;
    (*poos) << "ninternal: " << nintnodes_ << endl;
    (*poos) << "branch lengths: " << std::boolalpha << has_branchlengths_ << endl;
    if (has_branchlengths_) {
        if(rooted_tree_) {
            (*poos) << "rttipvar: " << rtvar_ << endl;
        } else {
            (*poos) << "rttipvar: NA" << endl;
        }
        (*poos) << "treelength: " << treelength_ << endl;
        (*poos) << "ultrametric: " << std::boolalpha << ultrametric_tree_ << endl;
    } else {
        (*poos) << "rttipvar: NA" << endl;
        (*poos) << "treelength: NA" << endl;
        (*poos) << "ultrametric: NA" << endl;
    }
    if (ultrametric_tree_) {
        (*poos) << "rootheight: " << rootheight_ << endl;
    } else {
        (*poos) << "rootheight: NA" << endl;
    }
}
