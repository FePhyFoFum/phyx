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
#include <numeric>
#include <algorithm>

using namespace std;

#include "node.h"
#include "tree.h"
#include "tree_utils.h"
#include "tree_reader.h"
#include "utils.h"

extern double EPSILON;


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


/*
 * calculate the variance between the lengths to the root and the tips
 */
double get_root_tip_var(Tree * tr){
    vector <double> paths;
    for (int i = 0; i < tr->getExternalNodeCount(); i++) {
        paths.push_back(get_length_to_root(tr->getExternalNode(i)));
    }
    double var = variance(paths);
    return var;
}


// assumes annotations are of form: [something]
void remove_annotations(Tree * tr) {
    for (int i=0; i < tr->getInternalNodeCount(); i++) {
        tr->getInternalNode(i)->setComment("");
    }
    for (int i=0; i < tr->getExternalNodeCount(); i++) {
        tr->getExternalNode(i)->setComment("");
    }
}


// same as above, but for names
void remove_internal_names(Tree * tr) {
    for (int i=0; i < tr->getInternalNodeCount(); i++) {
        tr->getInternalNode(i)->setName("");
    }
}


void remove_tips(Tree * tree, vector<string> & names, bool const& silent) {
    int num_names = names.size();
    int counter = 0;
    
    // new: note tree rooted status. if originally unrooted, make sure it stays that way on pruning
    bool rs = is_rooted(tree);
    //cout << "rooted tree: " << std::boolalpha << rs << endl;
    
    for (int i=0; i < num_names; i++) {
        //cout << "Attempting to remove tip '" << names[i] << "'." << endl;
        Node * m = tree->getExternalNode(names[i]);
        if (m != NULL) {
            tree->pruneExternalNode(m);
            counter++;
        } else {
            if (!silent) {
                cerr << names[i] << " not in tree"  << endl;
            }
        }
        //cout << "After pruning, tree rootedness: " << std::boolalpha << is_rooted(tree) << endl;
        if (rs != is_rooted(tree)) {
            //cout << "Unrooting tree..." << endl;
            // it is possible to go from unrooted to rooted on pruning, but not the other way (i think)
            tree->unRoot();
        }
        // debugging
        //cout << getNewickString(tree) << endl;
    }
}


void paint_nodes(Tree * tree, vector<string> & names, bool const& silent) {
    int num_names = names.size();
    tree->getRoot()->setPainted(true);
    for (int i=0; i < num_names; i++) {
        Node * m = tree->getNode(names[i]);
        if (m != NULL) {
            m->setPainted(true);
            Node * cur = m;
            while (cur != tree->getRoot()){
                cur = cur->getParent();
                if (cur->getPainted() == true)
                    break;
                else
                    cur->setPainted(true);
            }
        } else {
            if (!silent) {
                cerr << names[i] << " not in tree"  << endl;
            }
        }
    }
}


// treemap: key is focal node, value is vector of adjacent nodes
// both keys and values can be internal or terminal nodes
void create_tree_map_from_rootnode(Tree * tr, map<Node*,vector<Node*> > & tree_map) {
    
    bool debug = true;
    
    //check if rooted or unrooted
    bool rooted = is_rooted(tr);
    if (debug) cout << "tree is rooted: " << std::boolalpha << rooted << endl;
    for (int i=0; i < tr->getInternalNodeCount(); i++) {
        Node * tnd = tr->getInternalNode(i);
        if (debug) cout << "Focal node: " << tnd->getName() << endl;
        if (tnd->getParent() == NULL && rooted == true) { // root on rooted tree
            if (debug) cout << "  Node has no parent, as it is the root" << endl;
            continue;
        }
        vector<Node *> nds;
        for (int j=0; j < tnd->getChildCount(); j++) {
            nds.push_back(tnd->getChild(j));
            if (debug) cout << "  Adding child node: " << tnd->getChild(j)->getName() << endl;
        }
        if (tnd->getParent() == tr->getRoot() && rooted == true) {
            for (int j=0; j < tnd->getParent()->getChildCount(); j++) {
                if (tnd->getParent()->getChild(j) != tnd) {
                    nds.push_back(tnd->getParent()->getChild(j));
                    if (debug) cout << "  Adding sibling node: " << tnd->getChild(j)->getName() << endl;
                }
            }
        } else {
            if (tnd->getParent() != NULL) {
                nds.push_back(tnd->getParent());
                if (debug) cout << "  Adding parent node: " << tnd->getParent()->getName() << endl;
            }
        }
        tree_map[tnd] = nds;
    }
    for (int i=0; i < tr->getExternalNodeCount(); i++) {
        vector<Node *> nds;
        Node * tnd = tr->getExternalNode(i);
        if (debug) cout << "Focal node: " << tnd->getName() << endl;
        if (tnd->getParent() == tr->getRoot() && rooted == true) {
            for (int j=0; j < tnd->getParent()->getChildCount(); j++) {
                if (tnd->getParent()->getChild(j) != tnd) {
                    nds.push_back(tnd->getParent()->getChild(j));
                    if (debug) cout << "  Adding sibling node: " << tnd->getParent()->getChild(j)->getName() << endl;
                }
            }
        } else {
            nds.push_back(tnd->getParent());
            if (debug) cout << "  Adding parent node: " << tnd->getParent()->getName() << endl;
        }
        tree_map[tnd] = nds;
    }
    // print map<Node*, vector<Node*> >
    if (debug) {
        cout << endl << "TREE MAP:" << endl;
        for (map<Node*,vector<Node*> >::iterator it = tree_map.begin(); it != tree_map.end(); it++) {
            cout << "Node: " << it->first->getName() << endl;
            vector<Node*> terp = it->second;
            for (unsigned int i = 0; i < terp.size(); i++) {
                cout << "  " << terp[i]->getName() << " ";
            }
            cout << endl;
        }
    }
}


void nni_from_tree_map(Tree * tr, map<Node*,vector<Node*> > & tree_map) {
    
    bool debug = true;
    
    bool success = false;
    while (!success) {
        map<Node*,vector<Node*> >::iterator item = tree_map.begin();
        int r = random_int_range(0, tree_map.size());
        
        std::advance( item, r );
        Node * first = (*item).first;
        if (debug) cout << endl << "Node first (" << r << "): " << first->getName() << endl;
        
        // ack. 'middle' is not necessarily in the middle at all
        int r2 = random_int_range(0,(*item).second.size());
        item = tree_map.begin();
        std::advance( item, r2 );
        Node * middle = (*item).first;
        if (first == middle) {
            if (debug) cout << "Bailing because first == middle..." << endl;
            continue;
        }
        if (debug) cout << "Node middle (" << r2 << "): " << middle->getName() << endl;
        
        // furthermore, 'second' need not be anywhere near 'first' or 'middle
        int r3 = random_int_range(0,(*item).second.size());
        item = tree_map.begin();
        std::advance( item, r3 );
        Node * second = (*item).first;
        
        if (debug) cout << "Node second (" << r3 << "): " << second->getName() << endl;
        
        //TODO: need to fix what happens when the parent is the root, seems to break down
        if (first == second || second == middle || first == tr->getRoot()
            || second == tr->getRoot() || first->getParent() == tr->getRoot()
            || second->getParent() == tr->getRoot()) {
            if (debug) cout << "Bailing on this combination..." << endl;
            continue;
        }

        tr->exchangeNodes(first, second);
        tr->processRoot();

        success = true;
    }
    return;
}


// moving to own function, as is generally useful
// NOTE: this does not recognize trees where there is a polytomy at the root
bool is_rooted (Tree * tr) {
    bool rooted = false;
    if (tr->getRoot()->getChildCount() == 2) {
        rooted = true;
    }
    return rooted;
}


bool is_binary (Tree * tr) {
    bool binary = false;
    
    int nintnodes = tr->getInternalNodeCount();
    int ntips = tr->getExternalNodeCount();
    
    if (is_rooted(tr)) {
        if (nintnodes == (ntips - 1)) {
            binary = true;
        }
    } else {
        if (nintnodes == (ntips - 2)) {
            binary = true;
        }
    }
    return binary;
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


// NEED THIS?!? probably not. maybe doesn't hurt to have another interface?
bool has_branchlengths (Tree * tr) {
    bool gotem = tr->hasEdgeLengths();
    return gotem;
}


// simply multiply each edge by some scalar (determined somehow)
void rescale_tree(Tree * tr, double const& scalef) {
    int numNodes = tr->getNodeCount();
    for (int i = 0; i < numNodes; i++) {
        double terp = tr->getNode(i)->getBL();
        if (terp != 0.0) {
            double newlength = terp * scalef;
            tr->getNode(i)->setBL(newlength);
        }
    }
}


/* Two possible approaches:
 * 1) get get_length_to_root for each external node
 *    - easy, but inefficient
 * 2) do a postorder traversal
 *    - efficient, but harder
 *
*/
// need to check if BLs are present
bool is_ultrametric_paths (Tree * tr) {
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
    // this might be better, as it exits earlier on a fail.
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
    
    // alternate approach: look at variance of paths, compare to EPSILON
    double var = variance(paths);
    if (var < EPSILON) {
        ultrametric = true;
        // if ultrametric, go ahead and set node heights
        Node * root = tr->getRoot();
        set_node_heights(root);
    }
    
    //ultrametric = all_equal(paths);

    return ultrametric;
}


// ultrametricity needs to be checked before this
void set_node_heights (Node * node) {
    if (node == NULL) {
        return;
    }
    if (node->getChildCount() > 0) {
        vector <double> heights;
        //bool parentHeight = 0.0; // not used
        for (int i = 0; i < node->getChildCount(); i++) {
            set_node_heights(node->getChild(i));
            heights.push_back(node->getChild(i)->getBL() + node->getChild(i)->getHeight());
        }
        node->setHeight(heights[0]);
    }
    if (node->isExternal()) { // tip
        node->setHeight(0.0);
    }
}


// seems to suffer from rounding error (all_equal)
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
            //double parentHeight = 0.0; // not used
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


// single name version
bool check_name_against_tree(Tree * tr, string const& name) {
    bool allgood = true;
    Node * nd = tr->getExternalNode(name);
    if (nd == NULL) {
        allgood = false;
    }
    return allgood;
}


// if 'silent', don't complain if any outgroups are missing
bool reroot(Tree * tree, vector<string> const& outgroups, bool const& silent) {
    bool success = false;
    vector<string> outgr = outgroups;
    if (!silent) {
        if (!check_names_against_tree(tree, outgr)) {
            return false;
        }
    } else {
        outgr = get_names_in_tree(tree, outgr); // this is where things got deleted
        if (outgr.empty()) {
            return true;
        }
    }
    Node * m = tree->getMRCA(outgr);
    if (m == NULL) {
        return false;
    }
    if (m == tree->getRoot()) {
        // check to see if the outgroups are just the children of the root
        // if so, then do this
        // tree->rootWithRootTips(outgr); // this function does not exist
        // if not, then do this
        // tree->getInternalMRCA(outgr); // this function is only partly written
        
        // this can fail if there is a polytomy at the root which includes outgroups
        // if it fails, take complement of outgroup, root on that
        
        //cout << "trying to root on ingroup instead" << endl;
        vector <string> ingroup = get_complement_tip_set(tree, outgr);
        m = tree->getMRCA(ingroup);
        if (m == tree->getRoot()) {
            //cout << "doh: mrca of ingroup is the same node!" << endl;
            bool done = false;
            // root randomly, then do desired rooting
            Node * n = NULL;
            while (!done) {
                for (unsigned int i = 0; i < ingroup.size(); i++) {
                    n = tree->getExternalNode(ingroup[i]);
                    Node * p = n->getParent();
                    if (p != tree->getRoot()) {
                        //cout << "yay! can try with parent of " << ingroup[i] << endl;
                        done = true;
                        break;
                    } else {
                        //cout << "nope: can't use " << ingroup[i] << endl;
                    }
                    if (i == (ingroup.size() - 1)) {
                        done = true;
                    }
                }
            }
            success = tree->reRoot(n);
            string intermediate = getNewickString(tree);
            //cout << "intermediate result: " << intermediate << endl;
            //delete tree; // apparently necessary
            TreeReader tr;
            tree = tr.readTree(intermediate);
            m = tree->getMRCA(outgr);
        }
    }
    success = tree->reRoot(m);
    return success;
}


// return all names that are found in tree
vector <string> get_names_in_tree (Tree * tr, vector<string> const& names) {
    vector <string> matched;
    for (unsigned int i = 0; i < names.size(); i++) {
    //cout << "Checking name '" << names[i] << "'." << endl;
        Node * nd = tr->getExternalNode(names[i]);
        if (nd != NULL) {
            matched.push_back(names[i]);
        }
    }
    return matched;
}


vector <string> get_complement_tip_set (Tree * tr, vector<string> const& orig_names) {
    vector <string> comp = get_tip_labels(tr);
    std::vector<string>::iterator it;
    for (unsigned int i = 0; i < orig_names.size(); i++) {
        // remove bad names
        it = std::find(comp.begin(), comp.end(), orig_names[i]);
        if (it != comp.end()) {
            comp.erase(it);
        }
    }
    return comp;
}


// returns a sorted vector of all terminal labels
vector <string> get_tip_labels (Tree * tr) {
    vector <string> labels;
    for (int i = 0; i < tr->getExternalNodeCount(); i++) {
        labels.push_back(tr->getExternalNode(i)->getName());
    }
    sort(labels.begin(), labels.end());
    return labels;
}


// remove all knuckles in a tree
void deknuckle_tree (Tree * tree) {
    for (int i=0; i < tree->getInternalNodeCount(); i++) {
        Node * tnd = tree->getInternalNode(i);
        if (tnd->isKnuckle()) {
            //cout << tnd->getName() << " is a KNUCKLE!" << endl;
            // attempt to remove knuckle
            remove_knuckle(tnd);
        }
    }
}


// remove a individual knuckle
void remove_knuckle (Node * node) {
    if (node->isKnuckle()) {
        double el = node->getBL(); // edge length
        Node * cnd = node->getChild(0);
        el += cnd->getBL();
        Node * pnd = node->getParent();
        
        // delete node, set grandchild as child to grandparent, set el
        pnd->addChild(*cnd);
        cnd->setBL(el);
        pnd->removeChild(*node);
    }
}

// different version from directly accessing node/getNewick (which this does eventually)
// 1. check presence of edge lengths
// 2. includes semicolon in string
string getNewickString (Tree * tree) {
    bool bl = tree->hasEdgeLengths();
    string phy = tree->getRoot()->getNewick(bl) + ";";
    return phy;
}

// same as above but with objects
string getNewickString (Tree * tree, string obj) {
    bool bl = tree->hasEdgeLengths();
    string phy = tree->getRoot()->getNewick(bl,obj) + ";";
    return phy;
}

// not checking whether rooted in the first place, because registers as not
// e.g. ((((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0.5);
bool has_root_edge (Tree * tr) {
    bool rootEdge = false;
    int nchild = tr->getRoot()->getChildCount();
    //cout << "Root has " << nchild << " children." << endl;
    if (nchild == 1) {
        rootEdge = true;
    }
    return rootEdge;
}

string double_to_str(double d){
    size_t len = std::snprintf(0, 0, "%.16f", d);
    std::string s(len+1, 0);
    std::snprintf(&s[0], len+1, "%.16f", d);
    s.pop_back();
    s.erase(s.find_last_not_of('0') + 1, std::string::npos);
    if(s.back() == '.') {
        s.pop_back();
    }
    return s;
}