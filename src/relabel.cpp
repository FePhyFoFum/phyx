// rename taxon labels

#include <iostream>
#include <fstream>
#include <set>

using namespace std;

#include "relabel.h"
#include "tree.h"
#include "utils.h"


Relabel::Relabel (string & cnamesf, string nnamesf, bool const& verbose) {
    store_name_lists (cnamesf, nnamesf);
    verbose_ = verbose;
}

void Relabel::store_name_lists (string & cnamesf, string nnamesf) {
    vector <string> terp;
    string line;
    int ccount = 0;
    int ncount = 0;
    
    ifstream ifc(cnamesf.c_str());
    while (getline (ifc, line)) {
        if (!line.empty()) {
            ccount++;
            terp.push_back(line);
        }
    }
    ifc.close();
    old_names_ = terp;
    terp.clear();
    
    // check all current names are unique (otherwise won't work with existing map)
    set <string> orig(old_names_.begin(), old_names_.end());
    if (orig.size() < old_names_.size()) {
        cerr << "Oops! Current name list contains duplicates. Exiting." << endl;
        exit(0);
    }
    
    // TODO? clean names to make them jointly newick/nexus compliant
    // should we 'correct' invalid names? or leave it to user?
    // can use get_valid_newick_label or get_valid_nexus_label from utils.cpp
    // downside: if we change these, will not match what the user expects downstream
    // quotes should be ok if spaces are the issue
    
    ifstream ifn(nnamesf.c_str());
    while (getline (ifn, line)) {
        if (!line.empty()) {
            ncount++;
            terp.push_back(line);
        }
    }
    ifn.close();
    new_names_ = terp;
    
    if (ccount != ncount) {
        cerr << "Current (" << ccount << ") and new (" << ncount
            << ") lists differ in length. Exiting." << endl;
        exit(0);
    } else {
        num_tax_ = (int)old_names_.size();
        for (int i = 0; i < num_tax_; i++) {
            name_map_[old_names_[i]] = new_names_[i];
        }
    }
    /*
    for(map<string, string>::const_iterator it = name_map_.begin(); it != name_map_.end(); ++it) {
        std::cout << it->first << " " << it->second << endl;
    }
    */
}

// if verbose, will report the failed matches
void Relabel::relabel_tree (Tree * tr) {
    // keep track of matches
    set <string> orig(old_names_.begin(), old_names_.end());
    for (int i=0; i < tr->getExternalNodeCount(); i++) {
        string str = tr->getExternalNode(i)->getName();
        if (name_map_.find(str) != name_map_.end()) {
            //cout << "Tree label '" << str << "' found in name list!" << endl;
            tr->getExternalNode(i)->setName(name_map_[str]);
            orig.erase(str);
        } else {
            // see if it is quotes that is messing us up
            replace_all(str, "'", "");
            if (name_map_.find(str) != name_map_.end()) {
                //cout << "Found it (" << str << ") this time!" << endl;
                tr->getExternalNode(i)->setName(name_map_[str]);
                orig.erase(str);
            }
        }  
    }
    // do internal labels as well. be quieter here (don't expect internal nodes to have labels)
    for (int i=0; i < tr->getInternalNodeCount(); i++) {
        string str = tr->getInternalNode(i)->getName();
        if (str == "") {
            continue;
        }
        if (name_map_.find(str) != name_map_.end()) {
            //cout << "Tree label '" << str << "' found in name list!" << endl;
            tr->getInternalNode(i)->setName(name_map_[str]);
        } else {
            // see if it is quotes that is messing us up
            replace_all(str, "'", "");
            if (name_map_.find(str) != name_map_.end()) {
                //cout << "Found it this time!" << endl;
                tr->getInternalNode(i)->setName(name_map_[str]);
                orig.erase(str);
            }
        }  
    }
    
    // report failed matches (if present)
    if (orig.size() > 0) {
        if (verbose_) {
            cerr << "The following names to match were not found in the tree:" << endl;
            for (auto elem : orig) {
                cerr << elem << endl;
            }
        }
    }
}

// now returns a boolean indicating success
bool Relabel::relabel_sequence (Sequence & seq) {
    string str = seq.get_id();
    if (name_map_.find(str) != name_map_.end()) {
        seq.set_id(name_map_[str]);
        return true;
    } else {
        return false;
    }
}

set <string> Relabel::get_names_to_replace () {
    set <string> orig(old_names_.begin(), old_names_.end());
    return orig;
}
