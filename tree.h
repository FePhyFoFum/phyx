/*
 * tree.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef TREE_H_
#define TREE_H_

#include <string>
#include <vector>

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
	Node * getMRCATraverse(Node * curn1,Node * curn2);
	void setHeightFromRootToNode(Node & inNode, double newHeight);
	double getGreatestDistance(Node * inNode);

public:
	Tree();
	Tree(Node * root);

	void addExternalNode(Node * tn);
	void addInternalNode(Node * tn);
	void pruneExternalNode(Node * node);
	Node * getExternalNode(int num);
	Node * getExternalNode(string & name);
	Node * getInternalNode(int num);
	Node * getInternalNode(string & name);
	Node * getNode(int num);
	int getNodeCount();
	int getExternalNodeCount();
	int getInternalNodeCount();
	Node * getRoot();
	void setRoot(Node * inroot);
	void unRoot(Node & inroot);
	void reRoot(Node * inroot);
	void tritomyRoot(Node * toberoot);
	Node * getMRCA(vector<string> innodes);
	Node * getMRCA(vector<Node *> innodes);
	void processRoot();

	void setHeightFromRootToNodes();
	void setHeightFromTipToNodes();

	~Tree();
};

#endif /* TREE_H_ */
