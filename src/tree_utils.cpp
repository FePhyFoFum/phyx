/*
 * tree_utils.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <map>
#include <vector>
#include <iterator>
#include <functional>
#include <iostream>

using namespace std;

#include "node.h"
#include "tree.h"
#include "tree_utils.h"
#include "utils.h"

int get_distance_between_two_nodes(Tree * tr, Node * nd1, Node * nd2) {
    vector<Node *> vnd;
    vnd.push_back(nd1);
    vnd.push_back(nd2);
    Node * mrca = tr->getMRCA(vnd);
    int count = 0;
    Node * cur = nd1;
    while (cur != mrca) {
        cur = cur->getParent();
        count+= 1;
    }
    cur = nd2;
    while (cur != mrca) {
            cur = cur->getParent();
        count += 1;
    }
    return count;
}

/*
 * calculates the branch lengths to the root
 */
double get_length_to_root(Node * n) {
    double length = 0;
    while (n->getParent() != NULL) {
        length += n->getBL();
        n = n->getParent();
    }
    return length;
}

void create_tree_map_from_rootnode(Tree * tr, map<Node*,vector<Node*> > & tree_map) {
    //check if rooted or unrooted
    bool rooted = is_rooted(tr);
    for (int i=0; i < tr->getInternalNodeCount(); i++) {
        Node * tnd = tr->getInternalNode(i);
        if (tnd->getParent() == NULL && rooted == true) {
            continue;
        }
        vector<Node *> nds;
        for (int j=0;j<tnd->getChildCount();j++) {
            nds.push_back(tnd->getChild(j));
        }
        if (tnd->getParent() == tr->getRoot() && rooted == true) {
            for (int j=0; j < tnd->getParent()->getChildCount(); j++) {
                if (tnd->getParent()->getChild(j) != tnd) {
                    nds.push_back(tnd->getParent()->getChild(j));
                }
            }
        } else {
            if (tnd->getParent() != NULL) {
                nds.push_back(tnd->getParent());
            }
        }
        tree_map[tnd] = nds;
    }
    for (int i=0; i < tr->getExternalNodeCount(); i++) {
        vector<Node *> nds;
        Node * tnd = tr->getExternalNode(i);
        if (tnd->getParent() == tr->getRoot() && rooted == true) {
            for (int j=0; j < tnd->getParent()->getChildCount(); j++) {
                if (tnd->getParent()->getChild(j) != tnd) {
                    nds.push_back(tnd->getParent()->getChild(j));
                }
            }
        } else {
            nds.push_back(tnd->getParent());
        }
        tree_map[tnd] = nds;
    }    
}

void nni_from_tree_map(Tree * tr, map<Node*,vector<Node*> > & tree_map) {
    bool success = false;
    while (success == false) {
        map<Node*,vector<Node*> >::iterator item = tree_map.begin();
        int r = random_int_range(0,tree_map.size());
        std::advance( item, r );
        Node * first = (*item).first;

        int r2 = random_int_range(0,(*item).second.size());
        item = tree_map.begin();
        std::advance( item, r2 );
        Node * middle = (*item).first;
        if (first == middle) {
            continue;
        }

        int r3 = random_int_range(0,(*item).second.size());
        item = tree_map.begin();
        std::advance( item, r3);
        Node * second = (*item).first;
        //TODO: need to fix what happens when the parent is the root, seems to break down
        if (first == second || second == middle || first == tr->getRoot() || second == tr->getRoot()
           || first->getParent() == tr->getRoot() || second->getParent() == tr->getRoot()) {
            continue;
        }

        tr->exchangeNodes(first,second);
        tr->processRoot();

        success = true;
    }
    return;
}

// maybe move this elsewhere...
bool check_names_against_tree(Tree * tr, vector<string> names) {
    bool allgood = true;
    for (unsigned int i = 0; i < names.size(); i++) {
    //cout << "Checking name '" << names[i] << "'." << endl;
        Node * nd = tr->getExternalNode(names[i]);
        if (nd == NULL) {
            cout << "Taxon '" << names[i] << "' not found in tree." << endl;
            allgood = false;
        }
    }
    return allgood;
}

// moving to own function, as is generally useful
bool is_rooted (Tree * tr) {
    bool rooted = false;
    if (tr->getRoot()->getChildCount() == 2) {
        rooted = true;
    }
    return rooted;
}

// sum of edge lengths
double get_tree_length (Tree * tr) {
   double length = 0.0;
   //Node * root = tr->getRoot(); // not used
   int numNodes = tr->getNodeCount();
   for (int i = 0; i < numNodes; i++) {
       length += tr->getNode(i)->getBL();
   }
   return length;
}


/* Two possible approaches:
 * 1) get get_length_to_root for each external node
 *    - easy, but inefficient
 * 2) do a postorder traversal
 *    - efficient, but harder
 *
*/
// need to check if BLs are present
bool is_ultrametric_tips (Tree * tr) {
    bool ultrametric = false;
    if (!is_rooted(tr)) {
        return ultrametric;
    }
    if (get_tree_length(tr) == 0) {
        return ultrametric;
    }

    vector <double> paths;
    for (int i = 0; i < tr->getExternalNodeCount(); i++) {
        paths.push_back(get_length_to_root(tr->getExternalNode(i)));
        //cout << "Path: " << paths[i] << endl;
    }
    // compare against first
    /*
    if (all_of(paths.begin()+1, paths.end(), bind(essentially_equal, placeholders::_1, paths[0]))) {
        ultrametric = true;
    }
    */
    // this might be better, as it exists earlier on a fail.
    /*
    vector <double>::iterator it = find_if_not(paths.begin()+1, paths.end(), bind(essentially_equal, placeholders::_1, paths[0]));
    if (it == end(paths)) {
        ultrametric = true;
    }
    */

    /*
    for (unsigned int i = 1; i < paths.size(); i++) {
        cout << "Comparing " << paths[0] << " to " << paths[i] << " = " << (paths[0] - paths[i]) << "." << endl;
        essentially_equal(paths[0], paths[i]);
    }

    exit(0);
    */

    ultrametric = all_equal(paths);

    return ultrametric;
}

bool is_ultrametric_postorder (Tree * tr) {
    bool ultrametric = true;
    if (!is_rooted(tr)) {
        return ultrametric;
    }
    if (get_tree_length(tr) == 0) {
        return ultrametric;
    }
    Node * root = tr->getRoot();
    ultrametric = postorder_ultrametricity_check(root, ultrametric);

    return ultrametric;
}


// postorder
bool postorder_ultrametricity_check (Node * node, bool & ultrametric) {
    if (ultrametric) {
        if (node == NULL) {
            return ultrametric;
        }
        if (node->getChildCount() > 0) {
            vector <double> heights;
            //bool parentHeight = 0.0; // not used
            for (int i = 0; i < node->getChildCount(); i++) {
                postorder_ultrametricity_check(node->getChild(i), ultrametric);
                heights.push_back(node->getChild(i)->getBL() + node->getChild(i)->getHeight());
            }
            if (!all_equal(heights)) {
                ultrametric = false;
                return ultrametric;
            } else {
                node->setHeight(heights[0]);
            }
        }
        if (node->isExternal()) { // tip
            node->setHeight(0.0);
        }
    }
    return ultrametric;
}


bool reroot(Tree * tree, vector<string> & outgr) {
    if (!check_names_against_tree(tree, outgr)) {
        return false;
    }
    Node * m = tree->getMRCA(outgr);
    if (m == NULL) {
        return false;
    }
    if (m == tree->getRoot()) {
        //check to see if the outgroups are just the children of the root
        //if so, then do this
        //tree->rootWithRootTips(outgr);
        //if not, then do this
        tree->getInternalMRCA(outgr);
        return true;
    }
    bool success = tree->reRoot(m);
    return success;
}
