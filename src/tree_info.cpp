#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "tree_info.h"
#include "tree.h"
#include "node.h"
#include "tree_utils.h"


TreeInfo::TreeInfo (Tree * intree):rooted_tree_(false), ultrametric_tree_(false),
        binary_tree_(false), has_branchlengths_(false), treelength_(0.0),
        nintnodes_(0.0), ntips_(0.0), rootheight_(0.0), rtvar_(0.0) {
    tree_ = intree;
    calc_stats();
}


TreeInfo::TreeInfo (Tree * intree, const bool& ultracheck, const bool& binarycheck,
        const bool& agecheck, const bool& rootedcheck, const bool& ntipcheck,
        const bool& lengthcheck, const bool& namecheck, const bool& rtvarcheck,
        std::ostream* poos):nintnodes_(0.0) {
    tree_ = intree;
    if (ultracheck) {
        ultrametric_tree_ = is_ultrametric_paths(tree_);
        (*poos) << std::boolalpha << ultrametric_tree_ << std::endl;
    } else if (rtvarcheck) {
        has_branchlengths_ = tree_->hasEdgeLengths();
        if (has_branchlengths_) {
            rtvar_ = get_root_tip_var(tree_);
            (*poos) << rtvar_ << std::endl;
        } else {
            (*poos) << "NA" << std::endl;
        }
    } else if (binarycheck) {
        binary_tree_ = is_binary(tree_);
        (*poos) << std::boolalpha << binary_tree_ << std::endl;
    } else if (agecheck) {
        ultrametric_tree_ = is_ultrametric_paths(tree_);
        if (ultrametric_tree_) {
            rootheight_ = tree_->getRoot()->getHeight();
            (*poos) << rootheight_ << std::endl;
        } else {
            (*poos) << "NA" << std::endl;
        }
    } else if (rootedcheck) {
        rooted_tree_ = is_rooted(tree_);
        (*poos) << std::boolalpha << rooted_tree_ << std::endl;
    } else if (ntipcheck) {
        ntips_ = tree_->getExternalNodeCount();
        (*poos) << ntips_ << std::endl;
    } else if (lengthcheck) {
        has_branchlengths_ = tree_->hasEdgeLengths();
        if (has_branchlengths_) {
            treelength_ = get_tree_length(tree_);
            (*poos) << treelength_ << std::endl;
        } else {
            (*poos) << "NA" << std::endl;
        }
    } else if (namecheck) {
        tip_labels_ = get_tip_labels(tree_);
        for (const auto & tip_label : tip_labels_) {
            (*poos) << tip_label << std::endl;
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
    if (rooted_tree_) {
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


void TreeInfo::get_stats (std::ostream* poos) const {
    (*poos) << "rooted: " << std::boolalpha << rooted_tree_ << std::endl;
    (*poos) << "binary: " << std::boolalpha << binary_tree_ << std::endl;
    (*poos) << "nterminal: " << ntips_ << std::endl;
    (*poos) << "ninternal: " << nintnodes_ << std::endl;
    (*poos) << "branch lengths: " << std::boolalpha << has_branchlengths_ << std::endl;
    if (has_branchlengths_) {
        if (rooted_tree_) {
            (*poos) << "rttipvar: " << rtvar_ << std::endl;
        } else {
            (*poos) << "rttipvar: NA" << std::endl;
        }
        (*poos) << "treelength: " << treelength_ << std::endl;
        (*poos) << "ultrametric: " << std::boolalpha << ultrametric_tree_ << std::endl;
    } else {
        (*poos) << "rttipvar: NA" << std::endl;
        (*poos) << "treelength: NA" << std::endl;
        (*poos) << "ultrametric: NA" << std::endl;
    }
    if (ultrametric_tree_) {
        (*poos) << "rootheight: " << rootheight_ << std::endl;
    } else {
        (*poos) << "rootheight: NA" << std::endl;
    }
}
