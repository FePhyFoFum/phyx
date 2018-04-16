#ifndef _COLLAPSE_TREE_H_
#define _COLLAPSE_TREE_H_

using namespace std;

#include "tree.h"

/*
so the way i have this set up is for the class to only involve
setting up the collapsing scenario (e.g. threshold).
in other words, the tree is not a member of the class.
rather, a tree (or trees) can be thrown at it and processed.
this seems better to me, rather than multiple (identical)
class construction/destruction operations.
idunno. maybe this is not best.
*/
        
class Collapser {
private:
    float threshold_;
    bool scale_set_;
    
    void guess_scale (double & threshold);
    
public:
    Collapser (double & threshold);
    void collapse_edges (Tree * tr);
};

#endif /* _COLLAPSE_TREE_H_ */