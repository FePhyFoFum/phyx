#ifndef PX_POLYTOMY_TREE_H
#define PX_POLYTOMY_TREE_H

class Tree; // forward declaration

// do stuff with polytomies
        
class Polytomy {
private:

    std::vector<std::string> terminals_to_prune_;
    
public:
    explicit Polytomy (const long int& seed);
    void sample_polytomies (Tree * tr);
};

#endif /* PX_POLYTOMY_TREE_H */
