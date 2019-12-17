#ifndef _CL_TR_H_
#define _CL_TR_H_

class Tree; // forward declaration
class Node; // forward declaration

class CleanTree {
private:
    Tree* tree_;
    Node * root_;
    bool rooted_tree_;
    bool ultrametric_tree_;
    bool binary_tree_;
    bool has_branchlengths_;
    //double treelength_;
    void clean_all ();

public:
    CleanTree (Tree * intree);
};

#endif /* _CL_TR_H_ */
