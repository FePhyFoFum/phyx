/*
 * node.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef NODE_H_
#define NODE_H_

#include <string>
#include <vector>
#include <map>
using namespace std;


#include "node_object.h"
#include "vector_node_object.h"
#include "superdouble.h"

class Node{
private:
    double BL;//branch lengths
    double height; // could be from tip or from root
    int number;
    string name;
    Node * parent;
    vector<Node *> children;
    map<string,NodeObject *> assoc;
    string comment;

public:
    Node();
    Node(Node * parent);
    Node(double bl,int number,string name, Node * parent);


    int get_num_leaves();
    vector<Node*> get_leaves();

    vector<Node*> getChildren();
    bool isExternal();
    bool isInternal();
    bool isRoot();
    bool hasParent();
    void setParent(Node & p);
    int getNumber();
    void setNumber(int n);
    double getBL();
    void setBL(double bl);
    double getHeight();
    void setHeight(double he);
    bool hasChild(Node & test);
    bool addChild(Node & c);
    bool removeChild(Node & c);
    Node * getChild(int c);
    string getName();
    string getComment();
    void setName(string s);
    void setComment(string s);
    string getNewick(bool bl);
    string getNewick(bool bl,string obj);
    Node * getParent();
    int getChildCount();
    void assocObject(string name,NodeObject & obj);
    NodeObject * getObject(string name);


    VectorNodeObject<Superdouble> seg_sp_stoch_map_revB_time; //segment specific rev B, combining the tempA and the ENLT
    VectorNodeObject<Superdouble> seg_sp_stoch_map_revB_number; //segment specific rev B, combining the tempA and the ENLT
    ~Node();

};

#endif /* NODE_H_ */
