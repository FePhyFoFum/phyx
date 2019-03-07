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
double get_root_tip_var(Tree * tr);
void rescale_tree(Tree * tr, double const& scalef);

void remove_annotations(Tree * tr);
void remove_internal_names(Tree * tr);
void remove_tips(Tree * tree, vector<string> & names, bool const& silent);
void paint_nodes(Tree * tree, vector<string> & names, bool const& silent);

bool is_rooted (Tree * tr);
bool is_binary (Tree * tr);
double get_tree_length (Tree * tr);
bool has_branchlengths (Tree * tr);
bool is_ultrametric_paths (Tree * tr);
void set_node_heights (Node * node);
bool is_ultrametric_postorder (Tree * tr);
bool postorder_ultrametricity_check (Node * node, bool & ultrametric);
bool has_root_edge (Tree * tr);

bool check_names_against_tree(Tree * tr, vector<string> names);
bool check_name_against_tree(Tree * tr, string const& name);
bool reroot(Tree * tree, vector<string> const& outgroups, bool const& silent);
vector <string> get_names_in_tree(Tree * tr, vector<string> const& names);
vector <string> get_complement_tip_set (Tree * tr, vector<string> const& orig_names);
vector <string> get_tip_labels (Tree * tr);
void deknuckle_tree (Tree * tree);
void remove_knuckle (Node * node);
string getNewickString (Tree * tree);
string getNewickString (Tree * tree, string object);

string double_to_str(double d);
#endif /* _TREE_UTILS_H_ */
