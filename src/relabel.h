#ifndef _RELABEL_TREE_H_
#define _RELABEL_TREE_H_

#include <map>
#include <vector>
#include <set>

#include "tree.h"
#include "sequence.h"

class Relabel {
private:
    int num_tax_;
    std::vector<std::string> old_names_;
    std::vector<std::string> new_names_;
    std::map<std::string, std::string> name_map_;
    bool verbose_;
    
    void store_name_lists (std::string& cnamesf, std::string nnamesf);
    
public:
    Relabel (std::string& cnamesf, std::string nnamesf, const bool& verbose);
    void relabel_tree (Tree * tr);
    bool relabel_sequence (Sequence& seq);
    std::set<std::string> get_names_to_replace ();
};

#endif /* _RELABEL_TREE_H_ */
