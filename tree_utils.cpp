/*
 * tree_utils.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include "node.h"
#include "tree.h"
#include "tree_utils.h"

int get_distance_between_two_nodes(Tree * tr, Node * nd1, Node * nd2){
	vector<Node *> vnd;
	vnd.push_back(nd1);vnd.push_back(nd2);
	Node * mrca = tr->getMRCA(vnd);
	int count = 0;
	Node * cur = nd1;
	while (cur != mrca){
		cur = cur->getParent();
		count+= 1;
	}
	cur = nd2;
	while(cur != mrca){
		cur = cur->getParent();
		count += 1;
	}
	return count;
}

