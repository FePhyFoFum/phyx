#ifndef _LS_TR_H_
#define _LS_TR_H_

#include <string>

using namespace std;

#include "tree.h"

class TreeInfo {
private:
    Tree* tree_;
    bool rooted_tree_;
    bool ultrametric_tree_;
    bool binary_tree_;
    bool has_branchlengths_;
    double treelength_;
    double nintnodes_;
    double ntips_;
    double rootheight_;
    vector <string> tip_labels_;
    
    void calc_stats ();

public:
    TreeInfo (Tree * intree);
    TreeInfo (Tree * intree, bool const& ultracheck, bool const& binarycheck,
        bool const& agecheck, bool const& rootedcheck, bool const& ntipcheck,
        bool const& lengthcheck, bool const& namecheck, ostream* poos);
    void get_stats (ostream* poos);
};

#endif /* _LS_TR_H_ */
