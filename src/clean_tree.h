#ifndef PX_CL_TR_H
#define PX_CL_TR_H

class Tree; // forward declaration
class Node; // forward declaration

class CleanTree {
private:
    Tree* tree_;
    Node * root_;
    
    // the following are not currently being used
    //bool rooted_tree_;
    //bool ultrametric_tree_;
    //bool binary_tree_;
    //bool has_branchlengths_;
    bool remove_root_edge_;
    bool remove_labels_;
    bool remove_knuckles_;
    void clean_properties ();

public:
    CleanTree (Tree * intree, bool remove_root_edge, bool remove_labels,
            bool remove_knuckles);
};

#endif /* PX_CL_TR_H */
