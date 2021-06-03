#include <set>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <cassert>

#include "node.h"
#include "tree.h"


Tree::Tree ():internal_node_count_(0), external_node_count_(0), edge_lengths_(false),
        node_annotations_(false), internal_node_names_(false) {
    root_ = nullptr;
    processRoot();
}


Tree::Tree (Node * inroot):internal_node_count_(0), external_node_count_(0),
        edge_lengths_(false), node_annotations_(false), internal_node_names_(false) {
    root_ = inroot;
    processRoot();
}


void Tree::addExternalNode (Node * tn) {
    external_nodes_.push_back(tn);
    external_node_count_++;
    nodes_.push_back(tn);
}


void Tree::addInternalNode (Node * tn) {
    internal_nodes_.push_back(tn);
    internal_node_count_++;
    nodes_.push_back(tn);
}


Node * Tree::getExternalNode (int num) {
    return external_nodes_.at(static_cast<unsigned long>(num));
}


/*
 * could precompute this, check for run time differences
 */
Node * Tree::getExternalNode (const std::string& name) {
    Node * ret = nullptr;
    for (int i = 0; i < external_node_count_; i++) {
        if (external_nodes_.at(static_cast<unsigned long>(i))->getName() == name) {
            ret = external_nodes_.at(static_cast<unsigned long>(i));
        }
    }
    return ret;
}


Node * Tree::getInternalNode (int num) {
    return internal_nodes_.at(static_cast<unsigned long>(num));
}


/*
 * could precompute this, check for run time differences
 */
Node * Tree::getInternalNode (std::string& name) {
    Node * ret = nullptr;
    for (int i = 0; i < internal_node_count_; i++) {
        if (internal_nodes_.at(static_cast<unsigned long>(i))->getName() == name) {
            ret = internal_nodes_.at(static_cast<unsigned long>(i));
        }
    }
    return ret;
}


int Tree::getExternalNodeCount () const {
    return external_node_count_;
}


int Tree::getExtantNodeCount () {
    setHeightFromRootToNodes();
    double largest = 0;
    for (int i = 0; i < external_node_count_; i++) {
        double th = getExternalNode(i)->getHeight();
        if (th > largest) {
            largest = th;
        }
    }
    int count = 0;
    for (int i = 0; i < external_node_count_; i++) {
        double th = getExternalNode(i)->getHeight();
        if (fabs(th-largest) < 0.00001) {
            count += 1;
        }
    }
    return count;
}


// NOTE: this includes the root
int Tree::getInternalNodeCount () const {
    return internal_node_count_;
}


Node * Tree::getNode (int num) {
    return nodes_.at(static_cast<unsigned long>(num));
}


Node * Tree::getNode (std::string& name) {
    Node * ret = nullptr;
    if (name_node_map_.empty()) {
        for (auto & node : nodes_) {
            if (!node->getName().empty()) {
                name_node_map_[node->getName()] = node;
            }
        }
    }
    if (name_node_map_.count(name) != 0) {
        ret = name_node_map_[name];
    }
    return ret;
}


int Tree::getNodeCount () const {
    return static_cast<int>(nodes_.size());
}


Node * Tree::getRoot () const {
    return root_;
}


void Tree::setRoot (Node * inroot) {
    root_ = inroot;
}


void Tree::setEdgeLengthsPresent (bool res) {
    edge_lengths_ = res;
}


void Tree::setNodeAnnotationsPresent (bool res) {
    node_annotations_ = res;
}


bool Tree::hasNodeAnnotations () const {
    return node_annotations_;
}


void Tree::setNodeNamesPresent (bool res) {
    internal_node_names_ = res;
}


bool Tree::hasNodeNames () const {
    return internal_node_names_;
}


bool Tree::hasEdgeLengths () const {
    return edge_lengths_;
}


void Tree::unRoot () {
    processRoot();
    if (this->getRoot()->getChildCount() < 3) {
        tritomyRoot(nullptr);
        processRoot();
    }
}


/*
 * seems to be working but check for leaks
 */
bool Tree::reRoot (Node * inroot) {
    processRoot();
    bool ret = true;
    if (this->getRoot()->getChildCount() < 3) {
        tritomyRoot(nullptr); // not sure if this should actually be the inroot instead of nullptr
    }
    //std::cout << this->root->getNewick(false) << std::endl;
    if (root_ == inroot) {
        std::cerr << "you asked to root at the current root" << std::endl;
        ret = false;
    } else {
        Node * tempParent = inroot->getParent();
        auto * newRoot = new Node(tempParent);
        newRoot->addChild(*inroot);
        inroot->setParent(*newRoot);
        tempParent->removeChild(*inroot);
        tempParent->addChild(*newRoot);
        newRoot->setParent(*tempParent);
        newRoot->setBL(inroot->getBL() / 2);
        inroot->setBL(inroot->getBL() / 2);
        processReRoot(newRoot);
        setRoot(newRoot);
        processRoot();
        // duplicate support information (if present) for display purposes
        duplicateRootSupport();
    }
    return ret;
}


void Tree::removeRootEdge () {
    bool going = true; // multiple stems? nefarious...
    while (going) {
        if (root_->getChildCount() == 1) {
            /*
            1. grab child
            2. set child as new root
            3. delete old root
            4. reprocess tree?
            */
            Node * curr = root_->getChild(0);
            delete root_;
            setRoot(curr);
            processRoot();
        } else {
            break;
        }
    }
}


void Tree::duplicateRootSupport () {
    std::vector<Node*> kids = root_->getChildren();
    bool supfound = false;
    std::vector<std::string> sups;
    unsigned int numnodes = 0; // want to guard against when only 1 outgroup
    for (auto & kid : kids) {
        if (kid->isInternal()) {
            numnodes++;
            std::string x = kid->getName(); // support stored in name property
            if (!x.empty()) {
                supfound = true;
                sups.push_back(x);
            }
        }
    }
    if (supfound) {
        if (numnodes > sups.size()) {
            if (sups.size() == 1) {
                for (auto & kid : kids) {
                    if (kid->isInternal()) {
                        std::string x = kid->getName();
                        if (x.empty()) {
                            kid->setName(sups[0]);
                        }
                    }
                }
            } else {
                std::cout << "i don't know how this might happen..." << std::endl;
            }
        }
    }
}


/*
 * seems to be working now
 */
void Tree::tritomyRoot (Node * toberoot) {
    Node * curroot = this->getRoot();
    if (toberoot == nullptr) {
        if (curroot->getChild(0)->isInternal()) {
            Node * currootCH = curroot->getChild(0);
            double nbl = currootCH->getBL();
            curroot->getChild(1)->setBL(curroot->getChild(1)->getBL() + nbl);
            curroot->removeChild(*currootCH);
            for (int i = 0; i < currootCH->getChildCount(); i++) {
                curroot->addChild(*currootCH->getChild(i));
                //currootCH.getChild(i).setParent(curroot);
            }
        } else {
            Node * currootCH = curroot->getChild(1);
            double nbl = currootCH->getBL();
            curroot->getChild(0)->setBL(curroot->getChild(0)->getBL() + nbl);
            curroot->removeChild(*currootCH);
            for (int i = 0; i < currootCH->getChildCount(); i++) {
                curroot->addChild(*currootCH->getChild(i));
                //currootCH.getChild(i).setParent(curroot);
            }
        }
    } else {
        if (curroot->getChild(1) == toberoot) {
            Node * currootCH = curroot->getChild(0);
            double nbl = currootCH->getBL();
            curroot->getChild(1)->setBL(curroot->getChild(1)->getBL() + nbl);
            curroot->removeChild(*currootCH);
            for (int i = 0; i < currootCH->getChildCount(); i++) {
                curroot->addChild(*currootCH->getChild(i));
                //currootCH.getChild(i).setParent(curroot);
            }
        } else {
            Node * currootCH = curroot->getChild(1);
            double nbl = currootCH->getBL();
            curroot->getChild(0)->setBL(curroot->getChild(0)->getBL() + nbl);
            curroot->removeChild(*currootCH);
            for (int i = 0; i < currootCH->getChildCount(); i++) {
                curroot->addChild(*currootCH->getChild(i));
                //currootCH.getChild(i).setParent(curroot);
            }
        }
    }
}


Node * Tree::getMRCA (std::vector<std::string> innodes) {
    Node * mrca = nullptr;
    if (innodes.size() == 1) {
        return this->getExternalNode(innodes[0]);
    }
    std::vector<std::string> outgroup;
    outgroup.reserve(static_cast<unsigned int>(innodes.size()));
    for (const auto & innode : innodes) {
        outgroup.push_back(innode);
    }
    Node * cur1 = this->getExternalNode(outgroup.at(0));
    outgroup.erase(outgroup.begin());
    Node * cur2 = nullptr;
    Node * tempmrca = nullptr;
    while (!outgroup.empty()) {
        cur2 = this->getExternalNode(outgroup.at(0));
        outgroup.erase(outgroup.begin());
        tempmrca = getMRCATraverse(cur1, cur2);
        cur1 = tempmrca;
    }
    mrca = cur1;
    return mrca;
}


Node * Tree::getMRCA (std::vector<Node *> innodes) {
    Node * mrca = nullptr;
    if (innodes.size() == 1) {
        return innodes[0];
    }
    Node * cur1 = innodes.at(0);
    innodes.erase(innodes.begin());
    Node * cur2 = nullptr;
    Node * tempmrca = nullptr;
    while (!innodes.empty()) {
        cur2 = innodes.at(0);
        innodes.erase(innodes.begin());
        tempmrca = getMRCATraverse(cur1, cur2);
        cur1 = tempmrca;
    }
    mrca = cur1;
    return mrca;
}


/**
 * when the MRCA is returned as root, this will find 
 * the other node, internal that can serve as another root
 */
// looks like lots missing here!
Node * Tree::getInternalMRCA (std::vector<std::string>& innodes) {
    Node * mrca = nullptr;
    //std::set<Node *> original; // original set of nodes. not used
    if (innodes.size() == 1) {
        return this->getExternalNode(innodes[0]);
    } //else {
        // first get the root set
        // then get the leaves set
        // the difference will give the root back
        // we want to choose a node that has all of the outgroups
        // and the smallest number of other nodes
        // this is the best case when there
        for (int i = 0; i < this->getInternalNodeCount(); i++) {
            
            // NOTE: this is missing!
            
        }
    //}
    return mrca;
}


void Tree::setHeightFromRootToNodes () {
    setHeightFromRootToNode(*this->root_, this->root_->getBL());
}


// this should be called depth
void Tree::setHeightFromRootToNode (Node& inNode, double newHeight) {
    if (!inNode.isRoot()) {
        newHeight += inNode.getBL();
        inNode.setHeight(newHeight);
    } else {
        inNode.setHeight(newHeight);
    }
    for (int i = 0; i < inNode.getChildCount(); i++) {
        setHeightFromRootToNode(*inNode.getChild(i), newHeight);
    }
}


/*
 * only makes sense for ultrametric trees
 */
void Tree::setHeightFromTipToNodes () {
    for (int i = 0; i < external_node_count_; i++) {
        double curh = 0.0;
        Node * cur = this->getExternalNode(i);
        cur->setHeight(curh);
        while (cur->getParent() != nullptr) {
            curh += cur->getBL();
            cur = cur->getParent();
            if (cur->getHeight()<curh) {
                cur->setHeight(curh);
            }
        }
    }
}


/*
 * private
*/
void Tree::processRoot () {
    nodes_.clear();
    internal_nodes_.clear();
    external_nodes_.clear();
    internal_node_count_ = 0;
    external_node_count_ = 0;
    // this is always false
    /*
    if (&root_ == nullptr) {
        return;
    }
    */
    postOrderProcessRoot(root_);
}


void Tree::processReRoot (Node * node) {
    if (node->isRoot() || node->isExternal()) {
        return;
    }
    if (node->getParent() != nullptr) {
        processReRoot(node->getParent());
    }
    // Exchange branch label, length et cetera
    exchangeInfo(node->getParent(), node);
    // Rearrange topology
    Node * parent = node->getParent();
    node->addChild(*parent);
    parent->removeChild(*node);
    parent->setParent(*node);
}


void Tree::exchangeInfo(Node * node1, Node * node2) {
    std::string swaps;
    double swapd;
    swaps = node1->getName();
    node1->setName(node2->getName());
    node2->setName(swaps);

    swapd = node1->getBL();
    node1->setBL(node2->getBL());
    node2->setBL(swapd);
}


void Tree::exchangeNodes (Node * node1, Node * node2) {
    Node * par1 = node1->getParent();
    Node * par2 = node2->getParent();
    bool bp1 = par1->removeChild(*node1);
    assert(bp1);
    bool bp2 = par2->removeChild(*node2);
    assert(bp2);
    par1->addChild(*node2);
    par2->addChild(*node1);
}


void Tree::postOrderProcessRoot (Node * node) {
    if (node == nullptr) {
        return;
    }
    if (node->getChildCount() > 0) {
        for (int i = 0; i < node->getChildCount(); i++) {
            postOrderProcessRoot(node->getChild(i));
        }
    }
    if (node->isExternal()) {
        addExternalNode(node);
        node->setNumber(external_node_count_);
    } else {
        addInternalNode(node);
        node->setNumber(internal_node_count_);
    }
}


/*
    * how this works
    *
    * get the parent = parent
    * get the parent of the parent = mparent
    * remove parent from mparent
    * add !node from parent to mparent
    *
    * doesn't yet take care if node.parent == root
    * or polytomy
*/
void Tree::pruneExternalNode (Node * node) {
    if (node->isInternal()) {
        return;
    }
    Node * parent = node->getParent();
    if (parent->getChildCount() == 2) {
        Node * other = nullptr;
        for (int i = 0; i < parent->getChildCount(); i++) {
            if (parent->getChild(i) != node) {
                other = parent->getChild(i);
            }
        }
        double bl = other->getBL()+parent->getBL();
        Node * mparent = parent->getParent();
        if (mparent != nullptr) {
            mparent->addChild(*other);
            other->setBL(bl);
            for (int i = 0; i < mparent->getChildCount(); i++) {
                if (mparent->getChild(i) == parent) {
                    mparent->removeChild(*parent);
                    break;
                }
            }
        } else {
            root_ = other;
            //std::cout << "i am here apparently..." << std::endl;
        }
    } else {
        parent->removeChild(*node);
    }
    delete node;
    this->processRoot();
}


void Tree::pruneInternalNode (Node * node) {
    if (node->isExternal()) {
        return;
    }
    /*
     * how this works
     * 
     * store edge length of node: pel (points to `grandparent')
     * get the grandparent node = gparent
     * for each child of node, set parent to gparent, adding pel to existing edge length
     * delete node
     * 
     * doesn't yet deal with possibility that node is the root (i.e. no parent)
     * - this will arise with unrooted trees
    */
    double pel = node->getBL();
    //std::cout << "Parent edge length: " << pel << std::endl;
    Node * gparent = node->getParent();
    
    for (int i = 0; i < node->getChildCount(); i++) {
        Node * child = node->getChild(i);
        //std::cout << "  " << i << ". Original edge length: " << child->getBL() << std::endl;
        double newEL = pel + child->getBL();
        //std::cout << "     New edge length: " << newEL << std::endl;
        child->setBL(newEL);
        gparent->addChild(*child);
    }
    // now get rid of it
    gparent->removeChild(*node);
    delete node;
    this->processRoot();
}


Node * Tree::getMRCATraverse (Node * curn1, Node * curn2) {
    Node * mrca = nullptr;
    //get path to root for first node
    std::vector<Node *> path1;
    Node * parent = curn1;
    path1.push_back(parent);
    while (parent != nullptr) {
        path1.push_back(parent);
        if (parent->getParent() != nullptr) {
            parent = parent->getParent();
        } else {
            break;
        }
    }
    //find first match between this node and the first one
    parent = curn2;
    bool x = true;
    while (x) {
        for (auto & i : path1) {
            if (parent == i) {
                mrca = parent;
                x = false;
                break;
            }
        }
        parent = parent->getParent();
    }
    return mrca;
}


/*
 * end private
 */

Tree::~Tree () {
    for (int i = 0; i < internal_node_count_; i++) {
        delete getInternalNode(i);
    }
    for (int i = 0; i < external_node_count_; i++) {
        delete getExternalNode(i);
    }
}
