#include "tscale.h"
#include "tree.h"
#include "utils.h"
#include "tree_utils.h"


TScale::TScale ():scalef_(0.0), rootheight_(0.0), rootset_(false) {
}


void TScale::set_scalef (const double& scalef) {
    scalef_ = scalef;
    rootset_ = false;
}
 

void TScale::set_rootheight (const double& rootheight) {
    rootheight_ = rootheight;
    rootset_ = true;
}


void TScale::rescale (Tree * tr) {
    if (rootset_) {
        // need to figure out scaling factor from original and desired root ages
        // ultrametricity check is upstream
        double orig_rootheight = tr->getRoot()->getHeight();
        scalef_ = rootheight_ / orig_rootheight;
    }
    rescale_tree(tr, scalef_);
}
