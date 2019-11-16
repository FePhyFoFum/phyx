#ifndef _TSCALE_H_
#define _TSCALE_H_

#include "tree.h"

class TScale {
private:
    double scalef_;
    double rootheight_;
    bool rootset_;
    
public:
    TScale ();
    void set_scalef (const double& scalef);
    void set_rootheight (const double& rootheight);
    void rescale (Tree * tr);
};

#endif /* _TSCALE_H_ */
