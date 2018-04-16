// collapse internal edges below some support threshold

#include <iostream>
#include <math.h>       /* rint */

using namespace std;

#include "tree.h"
#include "tree_utils.h"
#include "collapse_tree.h"

Collapser::Collapser (double & threshold):scale_set_(false) {
    threshold_ = threshold;
    //cout << "Threshold set to: " << threshold_ << endl;
}


void Collapser::collapse_edges (Tree * tr) {
    //tree->hasNodeNames()
    
    bool isRooted = is_rooted(tr); // should be useful
    // cout << "isRooted: " << isRooted << endl;
    
    if (!tr->hasNodeAnnotations()) {
        //cout << "Dude. No annotations found in this tree. What are you even _doing_?!?" << endl;
    }
    if (!tr->hasNodeNames()) {
        //cout << "Dude. No node names found in this tree. What are you even _doing_?!?" << endl;
    } else {
        // if a node is removed, start over anew since root is reprocessed
        // should check if this is necessary, since could be expensive for large trees
        bool done = false;
        int loop_count = 0;
        while (!done) {
            loop_count++;
            //cout << "Loop #" << loop_count << ". There are " << tr->getInternalNodeCount() << " nodes to deal with." << endl;
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
                        break;
                    } else {
                        //cout << "This one is cool: " << cursup << endl;
                    }
                }
                if (i == (tr->getInternalNodeCount() - 1)) {
                    done = true;
                }
            }
        }
        //cout << "Went through loop " << loop_count << " times." << endl;
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

