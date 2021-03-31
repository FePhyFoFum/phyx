#ifndef _POLYTOMY_TREE_H_
#define _POLYTOMY_TREE_H_

class Tree; // forward declaration

// do stuff with polytomies
        
class Polytomy {
private:

    bool sample_polytomy_;
    std::vector<std::string> terminals_to_prune_;
    
public:
    Polytomy (const int& seed);
    void sample_polytomies (Tree * tr);
};

#endif /* _POLYTOMY_TREE_H_ */
