/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <deque>
#include <set>
#include <stack>

using namespace std;

#include "branch_segment.h"
#include "node.h"
#include "node_object.h"
#include "string_node_object.h"
#include "tree_utils.h"

Node::Node():BL(0.0), height(0.0), number(0), name(""), parent(NULL),
    children(vector<Node *> ()), assoc(map<string, NodeObject *>()),
    assocDV(map<string, vector<Superdouble> >()), comment(""), painted(false) {

}

Node::Node(Node * inparent):BL(0.0), height(0.0), number(0), name(""), parent(inparent),
    children(vector<Node *> ()), assoc(map<string,NodeObject *>()),
    assocDV(map<string, vector<Superdouble> >()), comment(""), painted(false) {

}

Node::Node(double bl, int innumber, string inname, Node * inparent):BL(bl), height(0.0),
    number(innumber), name(inname), parent(inparent), children(vector<Node *> ()),
    assoc(map<string, NodeObject *>()), assocDV(map<string,vector<Superdouble> >()), 
    comment(""), painted(false) {

}

vector<Node*> Node::getChildren() {
    return children;
}

bool Node::isExternal() {
    if (children.size() < 1) {
        return true;
    } else {
        return false;
    }
}

bool Node::isInternal() {
    if (children.size() > 0) {
        return true;
    } else {
        return false;
    }
}

bool Node::isRoot() {
    if (parent == NULL) {
        return true;
    } else {
        return false;
    }
}

// is the internal node of two degree (one parent and one child?)
bool Node::isKnuckle () {
    // neither tips nor the root can be a knuckle
    if (isRoot() || isExternal()) {
        return false;
    } else if (children.size() == 1) {
        return true;
    }
    return false;
}

bool Node::hasParent() {
    if (parent == NULL) {
        return false;
    } else {
        return true;
    }
}

void Node::setParent(Node & p) {
    parent = &p;
}

int Node::getNumber() {
    return number;
}

void Node::setNumber(int n) {
    number = n;
}

bool Node::getPainted(){ return painted;}
void Node::setPainted(bool p){ painted = p;}

double Node::getBL() {
    return BL;
}

void Node::setBL(double bl) {
    BL = bl;
}

double Node::getHeight() {
    return height;
}

void Node::setHeight(double he) {
    height = he;
}

double Node::getDepth() {
    return depth;
}

void Node::setDepth(double de) {
    depth = de;
}

bool Node::hasChild(Node & test) {
    bool ret = false;
    for (unsigned int i=0; i < children.size(); i++) {
        if (children.at(i) == &test) {
            ret = true;
            break;
        }
    }
    return ret;
}

bool Node::addChild(Node & c) {
    if (hasChild(c) == false) {
        children.push_back(&c);
        c.setParent(*this);
        return true;
    } else {
        return false;
    }
}

bool Node::removeChild(Node & c) {
    if (hasChild(c) == true) {
        for (unsigned int i=0; i < children.size(); i++) {
            if (children.at(i) == &c) {
                children.erase(children.begin()+i);
                break;
            }
        }
        return true;
    } else {
        return false;
    }
}

Node * Node::getChild(int c) {
    return children.at(c);
}

string Node::getName() {
    return name;
}

void Node::setName(string s) {
    name = s;
}

string Node::getComment() {
    return comment;
}

void Node::setComment(string s) {
    comment = s;
}

string Node::getNewick(bool bl) {
    string ret = "";
    for (int i=0; i < this->getChildCount(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret = ret+this->getChild(i)->getNewick(bl);
        if (bl == true) {
            //std::ostringstream o;
            ////20 is what you get from raxml
            //o.setf(ios::fixed,ios::floatfield);
            //o << setprecision(20) << this->getChild(i)->getBL();
            //ret += ":" + o.str();
            ret += ":" + double_to_str(this->getChild(i)->getBL());
        }
        if (i == this->getChildCount()-1) {
            ret += ")";
        } else {
            ret += ",";
        }
    }
    if (name.size() > 0) {
        ret += name;
    }
    return ret;
}

string Node::getPaintedNewick(bool bl){
    string ret = "";
    vector<int> paintedchildren;
    for (int i=0; i < this->getChildCount(); i++) {
        if (this->getChild(i)->getPainted() == true)
            paintedchildren.push_back(i);
    }
    for (unsigned int i=0; i < paintedchildren.size(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret += this->getChild(paintedchildren[i])->getPaintedNewick(bl);
        if (bl == true) {
            std::ostringstream o;
            //20 is what you get from raxml
            o << setprecision(20) << this->getChild(paintedchildren[i])->getBL();
            ret += ":" + o.str();
        }
        if (i == paintedchildren.size()-1) {
            ret += ")";
        } else {
            ret += ",";
        }
    }
    if (this->getName().size() > 0) {
        ret += this->getName();
    }
    return ret;
}

/*
 * should be returning the stringnodeobjects as the names for internal
 * nodes with the [&obj=string]
 * needs to be a string in setObject
 */
string Node::getNewick(bool bl, string obj) {
    string ret = "";
    for (int i=0; i < this->getChildCount(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret = ret+this->getChild(i)->getNewick(bl,obj);
        if (bl == true) {
            std::ostringstream o;
            o << this->getChild(i)->getBL();
            ret += ":" + o.str();
        }
        if (i == this->getChildCount()-1) {
            ret += ")";
        } else {
            ret += ",";
        }
    }
    if (this->name.size() > 0) {
        ret = ret + this->name;
    }
    if (this->getObject(obj) != NULL) {
        std::ostringstream o;
        o << (*((StringNodeObject*) (this->getObject(obj))));
        ret += "[&"+obj+"="+o.str()+"]";
    }
    return ret;
}

Node * Node::getParent() {
    return parent;
}

int Node::getChildCount() {
    return children.size();
}

void Node::assocObject(string name,NodeObject & obj) {
    if (assoc.count(name) > 0) {
        delete assoc[name];
    }
    assoc[name] = obj.clone();
}

void Node::assocDoubleVector(string name, vector<Superdouble> & obj) {
    if (assocDV.count(name) > 0) {
        assocDV.erase(name);
    }
    vector<Superdouble> tvec (obj.size());
    for (unsigned int i=0; i < obj.size(); i++) {
        tvec[i] = obj[i];
    }
    assocDV[name] = tvec;
}

vector<Superdouble> * Node::getDoubleVector(string name) {
    return &assocDV[name];
}

void Node::deleteDoubleVector(string name) {
    if (assocDV.count(name) > 0) {
        assocDV.erase(name);
    }
}

/*
 * gets the number of leaves from this node
 */
int Node::get_num_leaves() {
    int retnum = 0;
    stack<Node*> nodes;
    nodes.push(this);
    while (nodes.empty() == false) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i=0; i < nd->getChildCount(); i++) {
                nodes.push(nd->getChild(i));
            }
        } else {
            retnum += 1;
        }
    }
    return retnum;
}

/*
 * gets the leaves from this node
 */
vector<Node*> Node::get_leaves() {
    vector<Node*> retnodes;
    stack<Node*> nodes;
    nodes.push(this);
    while (nodes.empty() == false) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i=0; i < nd->getChildCount(); i++) {
                nodes.push(nd->getChild(i));
            }
        } else {
            retnodes.push_back(nd);
        }
    }
    return retnodes;
}

/*
 * gets the leaves from this node
 */
set<Node*> Node::get_leaves_set() {
    set<Node*> retnodes;
    stack<Node*> nodes;
    nodes.push(this);
    while (nodes.empty() == false) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i=0; i < nd->getChildCount(); i++) {
                nodes.push(nd->getChild(i));
            }
        } else {
            retnodes.insert(nd);
        }
    }
    return retnodes;
}

set<string> Node::get_leave_names_set() {
    stack<Node*> nodes;
    set<string> names;
    nodes.push(this);
    while (nodes.empty() == false) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i=0; i < nd->getChildCount(); i++) {
                nodes.push(nd->getChild(i));
            }
        } else {
            names.insert(nd->getName());
        }
    }
    return names;
}

vector<string> Node::get_leave_names() {
    stack<Node*> nodes;
    vector<string> names;
    nodes.push(this);
    while (nodes.empty() == false) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i=0; i < nd->getChildCount(); i++) {
                nodes.push(nd->getChild(i));
            }
        } else {
            names.push_back(nd->getName());
        }
    }
    return names;
}

/*
 * use the string ones like this
 * StringNodeObject sno("...a node object");
 * tree.getRoot()->assocObject("test",sno);
 * cout << *((StringNodeObject*) (tree.getRoot()->getObject("test"))) << endl;
 *
 * and the vector like
 * VectorNodeObject<int> vno;
 * vno.push_back(1);vno.push_back(2);
 * tree.getRoot()->assocObject("testvno",vno);
 * cout << ((VectorNodeObject<int> *) (tree.getRoot()->getObject("testvno")))->at(0) << endl;
 */

NodeObject  * Node::getObject(string name) {
    return assoc[name];
}

void Node::initSegVector() {
    segs = new vector<BranchSegment> ();
}

vector<BranchSegment> * Node::getSegVector() {
    return segs;
}

void Node::deleteSegVector() {
    delete segs;
}

/*
 * delete the node
 */
Node::~Node() {
    map<string,NodeObject *>::iterator it;
    for (it = assoc.begin() ; it != assoc.end(); it++) {
        delete assoc[it->first];
    }
}
