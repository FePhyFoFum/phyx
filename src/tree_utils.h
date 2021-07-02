#ifndef PX_TREE_UTILS_H
#define PX_TREE_UTILS_H

#include <string>
#include <vector>
#include <map>

class Tree; // forward declaration
class Node; // forward declaration

class Tree;

int get_distance_between_two_nodes (Tree * tr, Node * nd1, Node * nd2);
void create_tree_map_from_rootnode (Tree * tr, std::map<Node*, std::vector<Node*> >&);
void nni_from_tree_map (Tree * tr, std::map<Node*, std::vector<Node*> >&);
double get_length_to_root (Node * n);
double get_root_tip_var (Tree * tr);
void rescale_tree (Tree * tr, const double& scalef);

void remove_annotations (Tree * tr);
void remove_internal_names (Tree * tr);
void remove_tips (Tree * tree, std::vector<std::string>& names, const bool& silent);
Tree * get_induced_tree (Tree * tree, std::vector<std::string>& names, const bool& silent);
void paint_nodes (Tree * tree, std::vector<std::string>& names, const bool& silent);

bool is_monophyletic (Tree * tree, std::vector<std::string> names, const bool& skip_missing);
bool is_rooted (Tree * tr);
bool is_binary (Tree * tr);
double get_tree_length (Tree * tr);
bool has_branchlengths (Tree * tr);
bool is_ultrametric_paths (Tree * tr);
void set_node_heights (Node * node);
bool is_ultrametric_postorder (Tree * tr);
bool postorder_ultrametricity_check (Node * node, bool& ultrametric);
bool has_root_edge (Tree * tr);

bool check_names_against_tree (Tree * tr, std::vector<std::string> names);
bool check_name_against_tree (Tree * tr, const std::string& name);
bool reroot (Tree * tree, const std::vector<std::string>& outgroups, const bool& silent);
std::vector<std::string> get_names_in_tree (Tree * tr, const std::vector<std::string>& names);
std::vector<std::string> get_complement_tip_set (Tree * tr, const std::vector<std::string>& orig_names);
std::vector<std::string> get_names_in_tree_regex (Tree * tr, const std::string& pattern);
std::vector<std::string> get_tip_labels (Tree * tr);
void deknuckle_tree (Tree * tree);
void remove_knuckle (Node * node);
std::string getNewickString (Tree * tree);
std::string getNexusString (Tree * tree);
std::string getNewickString (Tree * tree, const std::string& obj);

std::string double_to_str (double d);

unsigned long int get_num_possible_trees (const unsigned int& n, const bool& rooted);

void get_terminal_children (Node * node, std::vector<Node *>& children);
void get_terminal_children (Node * node, std::vector<std::string>& children);
void get_all_descendant_nodes (Node * node, std::vector<Node *>& children);
#endif /* PX_TREE_UTILS_H */
