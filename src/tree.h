/*
 * tree.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef _TREE_H_
#define _TREE_H_

using namespace std;

#include "node.h"

class Tree{
private:
    Node * root;
    vector<Node *> nodes;
    vector<Node *> internalNodes;
    vector<Node *> externalNodes;
    int internalNodeCount;
    int externalNodeCount;
    
    void processReRoot(Node * node);
    void exchangeInfo(Node * node1, Node * node2);
    void postOrderProcessRoot(Node * node);
    Node * getMRCATraverse(Node * curn1, Node * curn2);
    void setHeightFromRootToNode(Node & inNode, double newHeight);
    double getGreatestDistance(Node * inNode);
    
public:
    Tree();
    Tree(Node * root);
    
    void addExternalNode(Node * tn);
    void addInternalNode(Node * tn);
    void pruneExternalNode(Node * node);
    Node * getExternalNode(int num);
    Node * getExternalNode(string name);
    Node * getInternalNode(int num);
    Node * getInternalNode(string & name);
    Node * getNode(int num);
    int getNodeCount();
    int getExtantNodeCount();
    int getExternalNodeCount();
    int getInternalNodeCount();
    Node * getRoot();
    void setRoot(Node * inroot);
    void unRoot();
    bool reRoot(Node * inroot);
    void duplicateRootSupport();
    void tritomyRoot(Node * toberoot);
    Node * getMRCA(vector<string> innodes);
    Node * getMRCA(vector<Node *> innodes);
    Node * getInternalMRCA(vector<string> & innodes);
    void processRoot();
    void exchangeNodes(Node * node1, Node * node2);
    
    void setHeightFromRootToNodes();
    void setHeightFromTipToNodes();
    
    ~Tree();
};

#endif /* _TREE_H_ */
