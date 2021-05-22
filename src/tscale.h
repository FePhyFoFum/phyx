#ifndef PX__TSCALE_H
#define PX__TSCALE_H

class Tree; // forward declaration

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

#endif /* PX__TSCALE_H */
