#include <string>
#include <iostream>
#include <random>

#include "tree.h"
#include "tree_utils.h"
#include "polytomy.h"
#include "utils.h"
#include "node.h"


Polytomy::Polytomy (const int& seed) {
    if (seed == -1) {
        srand(get_clock_seed());
    } else {
        srand(seed);
    }
}


// for each polytomy, remove all but 2 random descendants
void Polytomy::sample_polytomies (Tree * tr) {
    bool verbose = false;
    std::vector<std::string> tchildren;
    int numChildren = 0;
    
    // welp with multiple trees do not want to keep this
    terminals_to_prune_.clear();
    
    // currently works up through the tree
    // hrm maybe go from root down? would save on duplicates
    for (int i = 0; i < tr->getInternalNodeCount(); i++) {
        tchildren.clear();
        Node * m = tr->getInternalNode(i);
        std::string str = m->getName();
        
        numChildren = m->getChildCount();
        if (numChildren > 2) {
            // select random taxa to prune
            std::vector<int> terp = sample_without_replacement(numChildren, (numChildren - 2));
            if (verbose) {
                std::cout << "numChildren for node " << m->getName() << " = " << numChildren << std::endl;
                for (int p = 0; p < numChildren; p++) {
                    std::cout << p << ". " << m->getChild(p)->getName() << std::endl;
                }
            }
            
            for (unsigned int j = 0; j < terp.size(); j++) {
                Node * n = m->getChild(terp[j]);
                if (n->isExternal()) {
                    terminals_to_prune_.push_back(n->getName());
                    if (verbose) {
                        std::cout << "Random tip chosen to purge: " << n->getName() << std::endl;
                    }
                } else {
                    // if internal, grab all terminal descendants and add to purge list
                    get_terminal_children(n, tchildren);
                    if (verbose) {
                        std::cout << "welp dealing with an internal node to prune here..." << std::endl;
                    }
                    for (unsigned int k = 0; k < tchildren.size(); k++) {
                        terminals_to_prune_.push_back(tchildren[k]);
                        if (verbose) {
                            std::cout << "Descendant tip to be purged: " << tchildren[k] << std::endl;
                        }
                    }
                }
            }
        }
    }
    if (!terminals_to_prune_.empty()) {
        // get rid of any duplicates (can happen with nested polytomies)
        terminals_to_prune_ = get_unique_elements(terminals_to_prune_);
        if (verbose) {
            print_vector(terminals_to_prune_);
        }
        remove_tips(tr, terminals_to_prune_, true);
    }
}
