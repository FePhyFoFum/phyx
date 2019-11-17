#include <string>
#include <map>
#include <vector>
#include <iterator>
#include <functional>
#include <iostream>
#include <numeric>
#include <algorithm>

#include "node.h"
#include "tree.h"
#include "tree_utils.h"
#include "tree_reader.h"
#include "utils.h"

extern double EPSILON;


int get_distance_between_two_nodes (Tree * tr, Node * nd1, Node * nd2) {
    std::vector<Node *> vnd;
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
double get_length_to_root (Node * n) {
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
double get_root_tip_var (Tree * tr) {
    std::vector<double> paths;
    for (int i = 0; i < tr->getExternalNodeCount(); i++) {
        paths.push_back(get_length_to_root(tr->getExternalNode(i)));
    }
    double var = variance(paths);
    return var;
}


// assumes annotations are of form: [something]
void remove_annotations (Tree * tr) {
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


// this is _very_ slow compared to painting the induced tree (pxtrt)
void remove_tips (Tree * tree, std::vector<std::string> & names, const bool& silent) {
    int num_names = names.size();
    int counter = 0;
    
    // new: note tree rooted status. if originally unrooted, make sure it stays that way on pruning
    bool rs = is_rooted(tree);
    //std::cout << "rooted tree: " << std::boolalpha << rs << std::endl;
    
    for (int i=0; i < num_names; i++) {
        //std::cout << "Attempting to remove tip '" << names[i] << "'." << std::endl;
        Node * m = tree->getExternalNode(names[i]);
        if (m != NULL) {
            tree->pruneExternalNode(m);
            counter++;
        } else {
            if (!silent) {
                std::cerr << names[i] << " not in tree"  << std::endl;
            }
        }
        //std::cout << "After pruning, tree rootedness: " << std::boolalpha << is_rooted(tree) << std::endl;
        if (rs != is_rooted(tree)) {
            //std::cout << "Unrooting tree..." << std::endl;
            // it is possible to go from unrooted to rooted on pruning, but not the other way (i think)
            tree->unRoot();
        }
        //std::cout << getNewickString(tree) << std::endl;
    }
}


// tree must be rooted. assumes this is checked upstream (so alternative can be considered)
Tree * get_induced_tree (Tree * tree, std::vector<std::string> & names, const bool& silent) {
    paint_nodes(tree, names, silent);
    std::string tracetreenewick = tree->getRoot()->getPaintedNewick(true) + ";";
    Tree * indTree = read_tree_string(tracetreenewick);
    deknuckle_tree(indTree); // guaranteed to have knuckles atm
    return indTree;
}


// assumes a rooted tree
void paint_nodes (Tree * tree, std::vector<std::string> & names, const bool& silent) {
    int num_names = names.size();
    tree->getRoot()->setPainted(true); // probably do not want this, but mrca is expensive
    
    for (int i=0; i < num_names; i++) {
        Node * m = tree->getNode(names[i]);
        if (m != NULL) {
            m->setPainted(true);
            Node * cur = m;
            while (cur != tree->getRoot()) {
                cur = cur->getParent();
                // break early if encountering an existing painted path
                if (cur->getPainted() == true) {
                    break;
                } else {
                    cur->setPainted(true);
                }
            }
        } else {
            if (!silent) {
                std::cerr << names[i] << " not in tree"  << std::endl;
            }
        }
    }
}


// treemap: key is focal node, value is vectorof adjacent nodes
// both keys and values can be internal or terminal nodes
void create_tree_map_from_rootnode (Tree * tr, std::map<Node*,
        std::vector<Node*> > & tree_map) {
    
    bool debug = true;
    
    //check if rooted or unrooted
    bool rooted = is_rooted(tr);
    if (debug) std::cout << "tree is rooted: " << std::boolalpha << rooted << std::endl;
    for (int i=0; i < tr->getInternalNodeCount(); i++) {
        Node * tnd = tr->getInternalNode(i);
        if (debug) std::cout << "Focal node: " << tnd->getName() << std::endl;
        if (tnd->getParent() == NULL && rooted == true) { // root on rooted tree
            if (debug) std::cout << "  Node has no parent, as it is the root" << std::endl;
            continue;
        }
        std::vector<Node *> nds;
        for (int j=0; j < tnd->getChildCount(); j++) {
            nds.push_back(tnd->getChild(j));
            if (debug) std::cout << "  Adding child node: " << tnd->getChild(j)->getName() << std::endl;
        }
        if (tnd->getParent() == tr->getRoot() && rooted == true) {
            for (int j=0; j < tnd->getParent()->getChildCount(); j++) {
                if (tnd->getParent()->getChild(j) != tnd) {
                    nds.push_back(tnd->getParent()->getChild(j));
                    if (debug) std::cout << "  Adding sibling node: " << tnd->getChild(j)->getName() << std::endl;
                }
            }
        } else {
            if (tnd->getParent() != NULL) {
                nds.push_back(tnd->getParent());
                if (debug) std::cout << "  Adding parent node: " << tnd->getParent()->getName() << std::endl;
            }
        }
        tree_map[tnd] = nds;
    }
    for (int i=0; i < tr->getExternalNodeCount(); i++) {
        std::vector<Node *> nds;
        Node * tnd = tr->getExternalNode(i);
        if (debug) std::cout << "Focal node: " << tnd->getName() << std::endl;
        if (tnd->getParent() == tr->getRoot() && rooted == true) {
            for (int j=0; j < tnd->getParent()->getChildCount(); j++) {
                if (tnd->getParent()->getChild(j) != tnd) {
                    nds.push_back(tnd->getParent()->getChild(j));
                    if (debug) std::cout << "  Adding sibling node: " << tnd->getParent()->getChild(j)->getName() << std::endl;
                }
            }
        } else {
            nds.push_back(tnd->getParent());
            if (debug) std::cout << "  Adding parent node: " << tnd->getParent()->getName() << std::endl;
        }
        tree_map[tnd] = nds;
    }
    // print map<Node*, vector<Node*> >
    if (debug) {
        std::cout << std::endl << "TREE MAP:" << std::endl;
        for (std::map<Node*, std::vector<Node*> >::iterator it = tree_map.begin(); it != tree_map.end(); it++) {
            std::cout << "Node: " << it->first->getName() << std::endl;
            std::vector<Node*> terp = it->second;
            for (unsigned int i = 0; i < terp.size(); i++) {
                std::cout << "  " << terp[i]->getName() << " ";
            }
            std::cout << std::endl;
        }
    }
}


void nni_from_tree_map (Tree * tr, std::map<Node*, std::vector<Node*> > & tree_map) {
    bool debug = true;
    bool success = false;
    while (!success) {
        std::map<Node*, std::vector<Node*> >::iterator item = tree_map.begin();
        int r = random_int_range(0, tree_map.size());
        
        std::advance( item, r );
        Node * first = (*item).first;
        if (debug) std::cout << std::endl << "Node first (" << r << "): " << first->getName() << std::endl;
        
        // ack. 'middle' is not necessarily in the middle at all
        int r2 = random_int_range(0,(*item).second.size());
        item = tree_map.begin();
        std::advance( item, r2 );
        Node * middle = (*item).first;
        if (first == middle) {
            if (debug) std::cout << "Bailing because first == middle..." << std::endl;
            continue;
        }
        if (debug) std::cout << "Node middle (" << r2 << "): " << middle->getName() << std::endl;
        
        // furthermore, 'second' need not be anywhere near 'first' or 'middle
        int r3 = random_int_range(0,(*item).second.size());
        item = tree_map.begin();
        std::advance( item, r3 );
        Node * second = (*item).first;
        
        if (debug) std::cout << "Node second (" << r3 << "): " << second->getName() << std::endl;
        
        //TODO: need to fix what happens when the parent is the root, seems to break down
        if (first == second || second == middle || first == tr->getRoot()
            || second == tr->getRoot() || first->getParent() == tr->getRoot()
            || second->getParent() == tr->getRoot()) {
            if (debug) std::cout << "Bailing on this combination..." << std::endl;
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
void rescale_tree (Tree * tr, double const& scalef) {
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

    std::vector<double> paths;
    for (int i = 0; i < tr->getExternalNodeCount(); i++) {
        paths.push_back(get_length_to_root(tr->getExternalNode(i)));
        //std::cout << "Path: " << paths[i] << std::endl;
    }
    // compare against first
    /*
    if (all_of(paths.begin()+1, paths.end(), bind(essentially_equal, placeholders::_1, paths[0]))) {
        ultrametric = true;
    }
    */
    // this might be better, as it exits earlier on a fail.
    /*
    vector<double>::iterator it = find_if_not(paths.begin()+1, paths.end(), bind(essentially_equal, placeholders::_1, paths[0]));
    if (it == end(paths)) {
        ultrametric = true;
    }
    */

    /*
    for (unsigned int i = 1; i < paths.size(); i++) {
        std::cout << "Comparing " << paths[0] << " to " << paths[i] << " = " << (paths[0] - paths[i]) << "." << std::endl;
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
        std::vector<double> heights;
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
            std::vector<double> heights;
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


bool check_names_against_tree (Tree * tr, std::vector<std::string> names) {
    bool allgood = true;
    for (unsigned int i = 0; i < names.size(); i++) {
    //std::cout << "Checking name '" << names[i] << "'." << std::endl;
        Node * nd = tr->getExternalNode(names[i]);
        if (nd == NULL) {
            std::cout << "Taxon '" << names[i] << "' not found in tree." << std::endl;
            allgood = false;
        }
    }
    return allgood;
}


// single name version
bool check_name_against_tree (Tree * tr, const std::string& name) {
    bool allgood = true;
    Node * nd = tr->getExternalNode(name);
    if (nd == NULL) {
        allgood = false;
    }
    return allgood;
}


// if 'silent', don't complain if any outgroups are missing
bool reroot (Tree * tree, const std::vector<std::string>& outgroups, const bool& silent) {
    bool success = false;
    std::vector<std::string> outgr = outgroups;
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
        
        //std::cout << "trying to root on ingroup instead" << std::endl;
        std::vector<std::string> ingroup = get_complement_tip_set(tree, outgr);
        m = tree->getMRCA(ingroup);
        if (m == tree->getRoot()) {
            //std::cout << "doh: mrca of ingroup is the same node!" << std::endl;
            bool done = false;
            // root randomly, then do desired rooting
            Node * n = NULL;
            while (!done) {
                for (unsigned int i = 0; i < ingroup.size(); i++) {
                    n = tree->getExternalNode(ingroup[i]);
                    Node * p = n->getParent();
                    if (p != tree->getRoot()) {
                        //std::cout << "yay! can try with parent of " << ingroup[i] << std::endl;
                        done = true;
                        break;
                    } else {
                        //std::cout << "nope: can't use " << ingroup[i] << std::endl;
                    }
                    if (i == (ingroup.size() - 1)) {
                        done = true;
                    }
                }
            }
            success = tree->reRoot(n);
            std::string intermediate = getNewickString(tree);
            //std::cout << "intermediate result: " << intermediate << std::endl;
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
std::vector<std::string> get_names_in_tree (Tree * tr, std::vector<std::string> const& names) {
    std::vector<std::string> matched;
    for (unsigned int i = 0; i < names.size(); i++) {
    //std::cout << "Checking name '" << names[i] << "'." << std::endl;
        Node * nd = tr->getExternalNode(names[i]);
        if (nd != NULL) {
            matched.push_back(names[i]);
        }
    }
    return matched;
}


std::vector<std::string> get_complement_tip_set (Tree * tr, std::vector<std::string> const& orig_names) {
    std::vector<std::string> comp = get_tip_labels(tr);
    std::vector<std::string>::iterator it;
    for (unsigned int i = 0; i < orig_names.size(); i++) {
        // remove bad names
        it = std::find(comp.begin(), comp.end(), orig_names[i]);
        if (it != comp.end()) {
            comp.erase(it);
        }
    }
    return comp;
}


// returns a sorted vectorof all terminal labels
std::vector<std::string> get_tip_labels (Tree * tr) {
    std::vector<std::string> labels;
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
            //std::cout << tnd->getName() << " is a KNUCKLE!" << std::endl;
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
std::string getNewickString (Tree * tree) {
    bool bl = tree->hasEdgeLengths();
    std::string phy = tree->getRoot()->getNewick(bl) + ";";
    return phy;
}


// same as above but with objects
std::string getNewickString (Tree * tree, std::string obj) {
    bool bl = tree->hasEdgeLengths();
    std::string phy = tree->getRoot()->getNewick(bl,obj) + ";";
    return phy;
}


// not checking whether rooted in the first place, because registers as not
// e.g. ((((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0.5);
bool has_root_edge (Tree * tr) {
    bool rootEdge = false;
    int nchild = tr->getRoot()->getChildCount();
    //std::cout << "Root has " << nchild << " children." << std::endl;
    if (nchild == 1) {
        rootEdge = true;
    }
    return rootEdge;
}


std::string double_to_str (double d) {
    size_t len = snprintf(0, 0, "%.16f", d);
    std::string s(len+1, 0);
    snprintf(&s[0], len+1, "%.16f", d);
    s.pop_back();
    s.erase(s.find_last_not_of('0') + 1, std::string::npos);
    if (s.back() == '.') {
        s.pop_back();
    }
    return s;
}

unsigned long int get_num_possible_trees (int const& n, const bool& rooted) {
    return doublefactorial(2 * (n + rooted) - 5);
}
