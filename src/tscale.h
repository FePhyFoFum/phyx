#ifndef _TSCALE_H_
#define _TSCALE_H_

#include "tree.h"

using namespace std;

class TScale {
private:
    double scalef_;
    double rootheight_;
    bool rootset_;
    
public:
    TScale ();
    void set_scalef (double const& scalef);
    void set_rootheight (double const& rootheight);
    void rescale (Tree * tr);
};

#endif /* _TSCALE_H_ */