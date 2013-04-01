/*
 * tree_utils.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <map>
#include <vector>
#include <iterator>

#include "node.h"
#include "tree.h"
#include "tree_utils.h"
#include "utils.h"

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

/*
 * calculates the branch lengths to the root
 */
double get_length_to_root(Node * n){
	double length = 0;
	while (n->getParent() != NULL){
		length += n->getBL();
		n = n->getParent();
	}
	return length;
}

void create_tree_map_from_rootnode(Tree * tr, map<Node*,vector<Node*> > & tree_map){
    //check if rooted or unrooted
    bool rooted = false;
    if(tr->getRoot()->getChildCount() == 2){
	rooted = true;
    }
    for (int i=0;i<tr->getInternalNodeCount();i++){
	Node * tnd = tr->getInternalNode(i);
	if(tnd->getParent()==NULL && rooted == true){
	    continue;
	}
	vector<Node *> nds;
	for(int j=0;j<tnd->getChildCount();j++){
	    nds.push_back(tnd->getChild(j));
	}
	if(tnd->getParent() == tr->getRoot() && rooted == true){
	    for(int j=0;j<tnd->getParent()->getChildCount();j++){
		if(tnd->getParent()->getChild(j)!=tnd){
		    nds.push_back(tnd->getParent()->getChild(j));
		}
	    }
	}else{
	    if(tnd->getParent()!=NULL){
		nds.push_back(tnd->getParent());
	    }
	}
	tree_map[tnd] = nds;
    }
    for (int i=0;i<tr->getExternalNodeCount();i++){
	vector<Node *> nds;
	Node * tnd = tr->getExternalNode(i);
	if(tnd->getParent() == tr->getRoot() && rooted == true){
	    for(int j=0;j<tnd->getParent()->getChildCount();j++){
		if(tnd->getParent()->getChild(j)!=tnd){
		    nds.push_back(tnd->getParent()->getChild(j));
		}
	    }
	}else{
	    nds.push_back(tnd->getParent());
	}
	tree_map[tnd] = nds;
    }    
}

void nni_from_tree_map(Tree * tr, map<Node*,vector<Node*> > & tree_map){
    bool success = false;
    while(success == false){
	map<Node*,vector<Node*> >::iterator item = tree_map.begin();
	int r = random_int_range(0,tree_map.size());
	std::advance( item, r );
	Node * first = (*item).first;

	int r2 = random_int_range(0,(*item).second.size());
	item = tree_map.begin();
	std::advance( item, r2 );
	Node * middle = (*item).first;
	if (first == middle){
	    continue;
	}

	int r3 = random_int_range(0,(*item).second.size());
	item = tree_map.begin();
	std::advance( item, r3);
	Node * second = (*item).first;
	//TODO: need to fix what happens when the parent is the root, seems to break down
	if(first == second || second == middle || first == tr->getRoot() || second == tr->getRoot()
	   || first->getParent() == tr->getRoot() || second->getParent() == tr->getRoot()){
	    continue;
	}

	tr->exchangeNodes(first,second);
	tr->processRoot();

	success = true;
    }
    return;
}

// maybe move this elsewhere...
bool check_names_against_tree(Tree * tr, vector<string> names) {
	bool allgood = true;
	
	for (int i = 0; i < (int)names.size(); i++) {
	    //cout << "Checking name '" << names[i] << "'." << endl;
		Node * nd = tr->getExternalNode(names[i]);
		if (nd == NULL) {
			cout << "Taxon '" << names[i] << "' not found in tree." << endl;
			allgood = false;
		}
	}
	return allgood;
}
