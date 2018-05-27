#ifndef _RELABEL_TREE_H_
#define _RELABEL_TREE_H_

#include <map>
#include <vector>

using namespace std;

#include "tree.h"
#include "sequence.h"

class Relabel {
private:
    int num_tax_;
    vector <string> old_names_;
    vector <string> new_names_;
    map <string, string> name_map_;
    bool verbose_;
    
    void store_name_lists (string & cnamesf, string nnamesf);
    
public:
    Relabel (string & cnamesf, string nnamesf, bool const& verbose);
    void relabel_tree (Tree * tr);
    bool relabel_sequence (Sequence & seq);
    set <string> get_names_to_replace ();
};



#endif /* _RELABEL_TREE_H_ */