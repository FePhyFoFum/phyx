#ifndef _RELABEL_H_
#define _RELABEL_H_

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
    
    void store_name_lists (string & cnamesf, string nnamesf);
    
public:
    Relabel (string & cnamesf, string nnamesf);
    void relabel_tree (Tree * tr);
    void relabel_sequence (Sequence & seq);
};



#endif /* _RELABEL_H_ */