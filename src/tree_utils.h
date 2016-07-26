/*
 * tree_utils.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef _TREE_UTILS_H_
#define _TREE_UTILS_H_

using namespace std;

#include "node.h"

class Tree;

int get_distance_between_two_nodes(Tree * tr, Node * nd1, Node * nd2);
void create_tree_map_from_rootnode(Tree * tr, map<Node*, vector<Node*> > &);
void nni_from_tree_map(Tree * tr, map<Node*, vector<Node*> > &);
double get_length_to_root(Node * n);

void remove_annotations(Tree * tr);

bool is_rooted (Tree * tr);
bool is_binary (Tree * tr);
double get_tree_length (Tree * tr);
bool is_ultrametric_paths (Tree * tr);
void set_node_heights (Node * node);
bool is_ultrametric_postorder (Tree * tr);
bool postorder_ultrametricity_check (Node * node, bool & ultrametric);

bool check_names_against_tree(Tree * tr, vector<string> names);
bool reroot(Tree * tree, vector<string> & outgr, bool const& silent);
vector <string> get_names_in_tree(Tree * tr, vector<string> const& names);
vector <string> get_tip_labels (Tree * tr);

#endif /* _TREE_UTILS_H_ */
