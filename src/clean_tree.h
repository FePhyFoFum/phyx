#ifndef _CL_TR_H_
#define _CL_TR_H_

using namespace std;

#include "tree.h"

class CleanTree {
private:
    Tree* tree_;
    Node * root_;
    bool rooted_tree_;
    bool ultrametric_tree_;
    bool binary_tree_;
    bool has_branchlengths_;
    double treelength_;

    void clean_all ();

public:
    CleanTree (Tree * intree);
//    CleanTree (Tree * intree, bool const& ultracheck, bool const& binarycheck,
//        bool const& agecheck, bool const& rootedcheck, bool const& ntipcheck,
//        bool const& lengthcheck, bool const& namecheck, bool const& rtvarcheck,
//        ostream* poos);
//    void get_stats (ostream* poos);
};



#endif /* _CL_TR_H_ */
