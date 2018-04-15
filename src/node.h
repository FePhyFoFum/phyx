/*
 * node.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef _NODE_H_
#define _NODE_H_

#include <map>
#include <set>

using namespace std;

#include "branch_segment.h"
#include "vector_node_object.h"
#include "superdouble.h"

class Node{
private:
    double BL; // branch length, points to parent
    double height; // could be from tip or from root
    double depth;
    int number;
    string name;
    Node * parent;
    vector<Node *> children;
    map<string,NodeObject *> assoc;
    map<string, vector<Superdouble> > assocDV;
    vector<BranchSegment> * segs;
    string comment;
    bool painted;

public:
    Node();
    Node(Node * parent);
    Node(double bl, int number, string name, Node * parent);
    
    int get_num_leaves();
    vector<Node*> get_leaves();
    vector<string> get_leave_names();
    set<string> get_leave_names_set();
    set<Node *> get_leaves_set();
    
    vector<Node*> getChildren();
    bool isExternal();
    bool isInternal();
    bool isRoot();
    bool isKnuckle();
    bool hasParent();
    void setParent(Node & p);
    int getNumber();
    void setNumber(int n);
    bool getPainted();
    void setPainted(bool p);
    double getBL();
    void setBL(double bl);
    double getHeight();
    void setHeight(double he);
    double getDepth();
    void setDepth(double de);
    bool hasChild(Node & test);
    bool addChild(Node & c);
    bool removeChild(Node & c);
    Node * getChild(int c);
    string getName();
    string getComment();
    void setName(string s);
    void setComment(string s);
    string getNewick(bool bl);
    string getNewick(bool bl, string obj);
    string getPaintedNewick(bool bl);
    Node * getParent();
    int getChildCount();
    void assocObject(string name, NodeObject & obj);
    void assocDoubleVector(string name, vector<Superdouble> & obj);
    vector<Superdouble> * getDoubleVector(string name);
    void deleteDoubleVector(string name);
    NodeObject * getObject(string name);
    void initSegVector();
    vector<BranchSegment> * getSegVector();
    void deleteSegVector();
    
    VectorNodeObject<Superdouble> seg_sp_stoch_map_revB_time; //segment specific rev B, combining the tempA and the ENLT
    VectorNodeObject<Superdouble> seg_sp_stoch_map_revB_number; //segment specific rev B, combining the tempA and the ENLT
    ~Node();

};

#endif /* _NODE_H_ */
