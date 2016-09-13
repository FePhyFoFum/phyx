#ifndef _RENAME_H_
#define _RENAME_H_

#include <map>
#include <vector>

using namespace std;

#include "tree.h"

class Rename {
private:
    int num_tax_;
    vector <string> old_names_;
    vector <string> new_names_;
    map <string, string> name_map_;
    
    void store_name_lists (string & cnamesf, string nnamesf);
    
public:
    Rename (string & cnamesf, string nnamesf);
    void rename_tree (Tree * tr);
};














#endif /* _RENAME_H_ */