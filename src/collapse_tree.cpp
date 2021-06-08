#include <string>
#include <iostream>
#include <cmath>  /* rint */

#include "tree.h"
#include "tree_utils.h"
#include "collapse_tree.h"
#include "utils.h"
#include "node.h"


Collapser::Collapser (const double& threshold):scale_set_(false),has_labels_(false),
        has_annotations_(false) {
    threshold_ = threshold;
}


void Collapser::set_sup_string (const std::string& str) {
    sup_string_ = str;
}


/*
need to consider both node `names' _and_ `comments'
former can be newick or Nexus, latter are only Nexus (e.g. BEAST)
- annotations have certain strings identifying support:
    - 'posterior=', 'prob=', 'label='
*/
void Collapser::collapse_edges (Tree * tr) {
    has_labels_ = tr->hasNodeNames();
    has_annotations_ = tr->hasNodeAnnotations();
    //bool isRooted = is_rooted(tr); // should be useful
    // cout << "isRooted: " << isRooted << std::endl;
    if (!tr->hasNodeAnnotations()) {
        //std::cout << "Dude. No annotations found in this tree. What are you even _doing_?!?" << std::endl;
    }
    if (!tr->hasNodeNames()) {
        //std::cout << "Dude. No node names found in this tree. What are you even _doing_?!?" << std::endl;
    } else {
        // if a node is removed, start over anew since root is reprocessed
        // should check if this is necessary, since could be expensive for large trees
        bool done = false;
        int loop_count = 0;
        while (!done) {
            loop_count++;
            //std::cout << "Loop #" << loop_count << ". There are " << tr->getInternalNodeCount()
            //  << " internal nodes to consider." << std::endl;
            for (unsigned int i = 0; i < tr->getInternalNodeCount(); i++) {
                Node * m = tr->getInternalNode(i);
                std::string str = m->getName();
                if (str.empty()) {
                    //std::cout << "Whoops. This node has no support value." << std::endl;
                } else {
                    double cursup = std::stod(str);
                    if (!scale_set_) {
                        guess_scale(cursup);
                    }
                    if (cursup < threshold_) {
                        //std::cout << "Welp. This one has got to go (" << cursup << "); el = "
                        //  << m->getBL() << std::endl;
                        tr->pruneInternalNode(m);
                        break;
                    }
                    //std::cout << "This one is cool: " << cursup << "; el = " << m->getBL() << std::endl;
                }
                if (i == (tr->getInternalNodeCount() - 1)) {
                    done = true; // exit loop, no more polytomies remain
                }
            }
        }
        //std::cout << "Went through loop " << loop_count << " times." << std::endl;
    }
}


/*
need to consider ranges:
1) 0-100 (e.g., bootstraps)
2) 0.0-1.0 (probs)
- hrm: currently take the 'first' support value encountered.
  - could a tree have both support values above and below 1?
  - e.g. 74.3 and 0.9
  - let us assume: no!
  - therefore, no need for rintf
*/
// using a support value (the first encountered), determine whether scale is proportions or percentages
// if appropriate will reset threshold (e.g. from 0.5 to 50)
void Collapser::guess_scale (const float& sup) {
    if (sup > 1.0f) { // overkill
        //std::cout << "Ok, looks like a percentage here guys." << std::endl;
        threshold_ *= 100;
        //std::cout << "New threshold set to " << threshold_ << std::endl;
    } else {
        //std::cout << "Ok, looks like a proportion (probability?)." << std::endl;
    }
    scale_set_ = true;
}
