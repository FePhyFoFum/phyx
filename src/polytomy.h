#ifndef PX__POLYTOMY_TREE_H
#define PX__POLYTOMY_TREE_H

class Tree; // forward declaration

// do stuff with polytomies
        
class Polytomy {
private:

    bool sample_polytomy_;
    std::vector<std::string> terminals_to_prune_;
    
public:
    Polytomy (const long int& seed);
    void sample_polytomies (Tree * tr);
};

#endif /* PX__POLYTOMY_TREE_H */
