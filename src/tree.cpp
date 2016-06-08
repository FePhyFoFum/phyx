/*
 * tree.cpp
 *
 */

#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;

#include "node.h"
#include "tree.h"

Tree::Tree() {
    root = NULL;
    processRoot();
}

Tree::Tree(Node * inroot) {
    root = inroot;
    processRoot();
}

void Tree::addExternalNode(Node * tn) {
    externalNodes.push_back(tn);
    externalNodeCount++;
    nodes.push_back(tn);
}

void Tree::addInternalNode(Node * tn) {
    internalNodes.push_back(tn);
    internalNodeCount++;
    nodes.push_back(tn);
}

Node * Tree::getExternalNode(int num) {
    return externalNodes.at(num);
}

/*
 * could precompute this, check for run time differences
 */
Node * Tree::getExternalNode(string name) {
    Node * ret = NULL;
    for (int i=0; i < externalNodeCount; i++) {
        if (externalNodes.at(i)->getName() == name) {
            ret = externalNodes.at(i);
        }
    }
    return ret;
}

Node * Tree::getInternalNode(int num) {
    return internalNodes.at(num);
}

/*
 * could precompute this, check for run time differences
 */
Node * Tree::getInternalNode(string & name) {
    Node * ret = NULL;
    for (int i=0; i < internalNodeCount; i++) {
        if (internalNodes.at(i)->getName() == name) {
            ret = internalNodes.at(i);
        }
    }
    return ret;
}

int Tree::getExternalNodeCount() {
    return externalNodeCount;
}

int Tree::getExtantNodeCount() {
    setHeightFromRootToNodes();
    double largest = 0;
    for (int i=0; i < externalNodeCount; i++) {
        double th = getExternalNode(i)->getHeight();
        if (th > largest) {
            largest = th;
        }
    }
    int count = 0;
    for (int i=0; i < externalNodeCount; i++) {
        double th = getExternalNode(i)->getHeight();
        if (fabs(th-largest) < 0.00001) {
            count += 1;
        }
    }
    return count;
}

// NOTE: this includes the root
int Tree::getInternalNodeCount() {
    return internalNodeCount;
}

Node * Tree::getNode(int num) {
    return nodes.at(num);
}

int Tree::getNodeCount() {
    return nodes.size();
}

Node * Tree::getRoot() {
    return root;
}

void Tree::setRoot(Node * inroot) {
    root = inroot;
}

void Tree::unRoot() {
    processRoot();
    if (this->getRoot()->getChildCount() < 3) {
        tritomyRoot(NULL);
        processRoot();
    }
}

/*
 * seems to be working but check for leaks
 */
bool Tree::reRoot(Node * inroot) {
    processRoot();
    if (this->getRoot()->getChildCount() < 3) {
        tritomyRoot(NULL); // not sure if this should actually be the inroot instead of NULL
    }
    //cout << this->root->getNewick(false) << endl;
    if (root == inroot) {
        cerr << "you asked to root at the current root" << endl;
        return false;
    } else {
        Node * tempParent = inroot->getParent();
        Node * newRoot = new Node(tempParent);
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
        return true;
    }
}

/*
 * seems to be working now
 */
void Tree::tritomyRoot(Node * toberoot) {
    Node * curroot = this->getRoot();
    if (toberoot == NULL) {
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

Node * Tree::getMRCA(vector<string> innodes) {
    Node * mrca = NULL;
    if (innodes.size() == 1) {
        return this->getExternalNode(innodes[0]);
    } else {
        vector<string> outgroup;
        for (unsigned int i=0; i < innodes.size(); i++) {
            outgroup.push_back(innodes.at(i));
        }
        Node * cur1 = this->getExternalNode(outgroup.at(0));
        outgroup.erase(outgroup.begin());
        Node * cur2 = NULL;
        Node * tempmrca = NULL;
        while (outgroup.size() > 0) {
            cur2 = this->getExternalNode(outgroup.at(0));
            outgroup.erase(outgroup.begin());
            tempmrca = getMRCATraverse(cur1, cur2);
            cur1 = tempmrca;
        }
        mrca = cur1;
    }
    return mrca;
}

Node * Tree::getMRCA(vector<Node *> innodes) {
    Node * mrca = NULL;
    if (innodes.size() == 1) {
        return innodes[0];
    } else {
        Node * cur1 = innodes.at(0);
        innodes.erase(innodes.begin());
        Node * cur2 = NULL;
        Node * tempmrca = NULL;
        while (innodes.size() > 0) {
            cur2 = innodes.at(0);
            innodes.erase(innodes.begin());
            tempmrca = getMRCATraverse(cur1, cur2);
            cur1 = tempmrca;
        }
        mrca = cur1;
    }
    return mrca;
}

/**
 * when the MRCA is returned as root, this will find 
 * the other node, internal that can serve as another root
 */
Node * Tree::getInternalMRCA(vector<string> & innodes) {
    Node * mrca = NULL;
    set<Node *> original;//original set of nodes
    if (innodes.size() == 1) {
        return this->getExternalNode(innodes[0]);
    } else {
        //first get the root set
        //then get the leaves set
        //the difference will give the root back
        //we want to choose a node that has all of the outgroups
        //and the smallest number of other nodes
        //this is the best case when there
        for (int i=0; i < this->getInternalNodeCount(); i++) {

        }
    }
    return mrca;
}

void Tree::setHeightFromRootToNodes() {
    setHeightFromRootToNode(*this->root,this->root->getBL());
}

void Tree::setHeightFromRootToNode(Node & inNode, double newHeight) {
    if (inNode.isRoot() == false) {
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
void Tree::setHeightFromTipToNodes() {
    for (int i = 0; i < externalNodeCount; i++) {
        double curh = 0.0;
        Node * cur = this->getExternalNode(i);
        cur->setHeight(curh);
        while (cur->getParent() != NULL) {
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
void Tree::processRoot() {
    nodes.clear();
    internalNodes.clear();
    externalNodes.clear();
    internalNodeCount = 0;
    externalNodeCount = 0;
    if (&root == NULL) {
        return;
    }
    postOrderProcessRoot(root);
}

void Tree::processReRoot(Node * node) {
    if (node->isRoot() || node->isExternal()) {
        return;
    }
    if (node->getParent() != NULL) {
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
    string swaps;
    double swapd;
    swaps = node1->getName();
    node1->setName(node2->getName());
    node2->setName(swaps);

    swapd = node1->getBL();
    node1->setBL(node2->getBL());
    node2->setBL(swapd);
}

void Tree::exchangeNodes(Node * node1, Node * node2) {
    Node * par1 = node1->getParent();
    Node * par2 = node2->getParent();
    bool bp1 = par1->removeChild(*node1);
    assert(bp1);
    bool bp2 = par2->removeChild(*node2);
    assert(bp2);
    par1->addChild(*node2);
    par2->addChild(*node1);
}

void Tree::postOrderProcessRoot(Node * node) {
    if (node == NULL)
        return;
    if (node->getChildCount() > 0) {
        for (int i = 0; i < node->getChildCount(); i++) {
            postOrderProcessRoot(node->getChild(i));
        }
    }
    if (node->isExternal()) {
        addExternalNode(node);
        node->setNumber(externalNodeCount);
    } else {
        addInternalNode(node);
        node->setNumber(internalNodeCount);
    }
}

void Tree::pruneExternalNode(Node * node) {
    if (node->isInternal()) {
        return;
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
    double bl = 0;
    Node * parent = node->getParent();
    Node * other = NULL;
    for (int i=0; i < parent->getChildCount(); i++) {
        if (parent->getChild(i) != node) {
            other = parent->getChild(i);
        }
    }
    bl = other->getBL()+parent->getBL();
    Node * mparent = parent->getParent();
    if (mparent != NULL) {
        mparent->addChild(*other);
        other->setBL(bl);
        for (int i=0; i < mparent->getChildCount(); i++) {
            if (mparent->getChild(i)==parent) {
                mparent->removeChild(*parent);
                break;
            }
        }
    } else {
        root = other;
    }
    delete node;
    this->processRoot();
}

Node * Tree::getMRCATraverse(Node * curn1, Node * curn2) {
    Node * mrca = NULL;
    //get path to root for first node
    vector<Node *> path1;
    Node * parent = curn1;
    path1.push_back(parent);
    while (parent != NULL) {
        path1.push_back(parent);
        if (parent->getParent() != NULL) {
            parent = parent->getParent();
        } else {
            break;
        }
    }
    //find first match between this node and the first one
    parent = curn2;
    bool x = true;
    while (x == true) {
        for (unsigned int i = 0; i < path1.size(); i++) {
            if (parent == path1.at(i)) {
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

Tree::~Tree() {
    for (int i=0; i < internalNodeCount; i++) {
        delete getInternalNode(i);
    }
    for (int i=0; i < externalNodeCount; i++) {
        delete getExternalNode(i);
    }
}
