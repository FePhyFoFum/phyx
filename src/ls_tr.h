#ifndef _LS_TR_H_
#define _LS_TR_H_

#include <string>

using namespace std;

#include "tree.h"

class TreeInfo {
private:
    string model;
    Tree* tree;
    bool rooted_tree;
    bool ultrametric_tree;
    bool binary_tree;
    bool has_branchlengths;
    double treelength;
    double nintnodes;
    double ntips;
    double rootheight;
    
    void calc_stats ();

public:
    TreeInfo (Tree * intree);
    TreeInfo (Tree * intree, bool const& ultracheck, bool const& binarycheck,
        bool const& agecheck, bool const& rootedcheck, bool const& ntipcheck,
        bool const& lengthcheck, ostream* poos);
    void get_stats (ostream* poos);
};

#endif /* _LS_TR_H_ */
