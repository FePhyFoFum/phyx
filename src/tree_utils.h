/*
 * tree_utils.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef _TREE_UTILS_H_
#define _TREE_UTILS_H_


#include "node.h"

class Tree;

int get_distance_between_two_nodes(Tree * tr, Node * nd1, Node * nd2);
void create_tree_map_from_rootnode(Tree * tr, std::map<Node*, std::vector<Node*> > &);
void nni_from_tree_map(Tree * tr, std::map<Node*, std::vector<Node*> > &);
double get_length_to_root(Node * n);
double get_root_tip_var(Tree * tr);
void rescale_tree(Tree * tr, double const& scalef);

void remove_annotations(Tree * tr);
void remove_internal_names(Tree * tr);
void remove_tips(Tree * tree, std::vector<std::string> & names, bool const& silent);
Tree * get_induced_tree (Tree * tree, std::vector<std::string> & names, bool const& silent);
void paint_nodes(Tree * tree, std::vector<std::string> & names, bool const& silent);

bool is_rooted (Tree * tr);
bool is_binary (Tree * tr);
double get_tree_length (Tree * tr);
bool has_branchlengths (Tree * tr);
bool is_ultrametric_paths (Tree * tr);
void set_node_heights (Node * node);
bool is_ultrametric_postorder (Tree * tr);
bool postorder_ultrametricity_check (Node * node, bool & ultrametric);
bool has_root_edge (Tree * tr);

bool check_names_against_tree(Tree * tr, std::vector<std::string> names);
bool check_name_against_tree(Tree * tr, std::string const& name);
bool reroot(Tree * tree, std::vector<std::string> const& outgroups, bool const& silent);
std::vector <std::string> get_names_in_tree(Tree * tr, std::vector<std::string> const& names);
std::vector <std::string> get_complement_tip_set (Tree * tr, std::vector<std::string> const& orig_names);
std::vector <std::string> get_tip_labels (Tree * tr);
void deknuckle_tree (Tree * tree);
void remove_knuckle (Node * node);
std::string getNewickString (Tree * tree);
std::string getNewickString (Tree * tree, std::string object);

std::string double_to_str(double d);

unsigned long int get_num_possible_trees (int const& n, bool const& rooted);
#endif /* _TREE_UTILS_H_ */
