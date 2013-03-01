/*
 * tree_utils.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef TREE_UTILS_H_
#define TREE_UTILS_H_

#include <map>
#include <vector>

#include "node.h"
#include "tree.h"

int get_distance_between_two_nodes(Tree * tr,Node * nd1, Node * nd2);
void create_tree_map_from_rootnode(Tree * tr, map<Node*,vector<Node*> > &);
void nni_from_tree_map(Tree * tr, map<Node*,vector<Node*> > &);

#endif /* TREE_UTILS_H_ */
