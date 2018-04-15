// collapse internal edges below some support threshold

#include <iostream>
#include <math.h>       /* rint */

using namespace std;

#include "tree.h"
#include "utils.h"
#include "collapse_tree.h"

Collapser::Collapser (double & threshold) {
    threshold_ = threshold;
    //cout << "Threshold set to: " << threshold_ << endl;
}


void Collapser::collapse_edges (Tree * tr) {
    //tree->hasNodeNames()
    if (!tr->hasNodeAnnotations()) {
        //cout << "Dude. No annotations found in this tree. What are you even _doing_?!?" << endl;
    }
    if (!tr->hasNodeNames()) {
        //cout << "Dude. No node names found in this tree. What are you even _doing_?!?" << endl;
    } else {
        for (int i=0; i < tr->getInternalNodeCount(); i++) {
            Node * m = tr->getInternalNode(i);
            string str = m->getName();
            if (str == "") {
                //cout << "Whoops. This node has no support value." << endl;
            } else {
                float cursup = stof(str);
                if (cursup < threshold_) {
                    //cout << "Welp. This one has got to go (" << cursup << ")." << endl;
                    tr->pruneInternalNode(m);
                } else {
                    //cout << "This one is cool: " << cursup << endl;
                }
            }
        }
    }
}



/*
need to consider both node `names' _and_ `comments'
former can be newick or Nexus, latter are only Nexus (e.g. BEAST)
*/


/*
need to consider ranges:
1) 0-100 (e.g., bootstraps)
2) 0.0-1.0 (probs)
- use rint to find out if what is passed in is an integer or float:
- rintf(f)==f
*/

void Collapser::guess_scale (double & threshold) {
    
}

