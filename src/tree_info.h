#ifndef PX_LS_TR_H
#define PX_LS_TR_H

#include <string>
#include <vector>
#include <iostream>

class Tree; // forward declaration

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
    double rtvar_;
    std::vector<std::string> tip_labels_;
    
    void calc_stats ();

public:
    TreeInfo (Tree * intree);
    TreeInfo (Tree * intree, const bool& ultracheck, const bool& binarycheck,
        const bool& agecheck, const bool& rootedcheck, const bool& ntipcheck,
        const bool& lengthcheck, const bool& namecheck, const bool& rtvarcheck,
        std::ostream* poos);
    void get_stats (std::ostream* poos);
};

#endif /* PX_LS_TR_H */
