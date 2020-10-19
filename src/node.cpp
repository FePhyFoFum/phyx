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

#include "branch_segment.h"
#include "node.h"
#include "node_object.h"
#include "string_node_object.h"
#include "tree_utils.h"
#include "utils.h"


Node::Node ():BL(0.0), height(0.0), number(0), name(""), parent(NULL),
    children(std::vector<Node *> ()), assoc(std::map<std::string, NodeObject *>()),
    assocDV(std::map<std::string, std::vector<Superdouble> >()), comment(""), painted(false) {

}


Node::Node (Node * inparent):BL(0.0), height(0.0), number(0), name(""), parent(inparent),
    children(std::vector<Node *> ()), assoc(std::map<std::string, NodeObject *>()),
    assocDV(std::map<std::string, std::vector<Superdouble> >()), comment(""), painted(false) {

}


Node::Node (double bl, int innumber, std::string inname, Node * inparent):BL(bl), height(0.0),
    number(innumber), name(inname), parent(inparent), children(std::vector<Node *> ()),
    assoc(std::map<std::string, NodeObject *>()), assocDV(std::map<std::string, std::vector<Superdouble> >()), 
    comment(""), painted(false) {

}


std::vector<Node*> Node::getChildren () {
    return children;
}


bool Node::isExternal () {
    if (children.size() < 1) {
        return true;
    } else {
        return false;
    }
}


bool Node::isInternal () {
    if (children.size() > 0) {
        return true;
    } else {
        return false;
    }
}


bool Node::isRoot () {
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


bool Node::hasParent () {
    if (parent == NULL) {
        return false;
    } else {
        return true;
    }
}


void Node::setParent (Node& p) {
    parent = &p;
}


int Node::getNumber () {
    return number;
}


void Node::setNumber (int n) {
    number = n;
}


bool Node::getPainted () {
    return painted;
}


void Node::setPainted (bool p) {
    painted = p;
}


double Node::getBL () {
    return BL;
}


void Node::setBL (double bl) {
    BL = bl;
}


double Node::getHeight () {
    return height;
}


void Node::setHeight (double he) {
    height = he;
}


double Node::getDepth () {
    return depth;
}


void Node::setDepth (double de) {
    depth = de;
}


bool Node::hasChild (Node& test) {
    bool ret = false;
    for (unsigned int i=0; i < children.size(); i++) {
        if (children.at(i) == &test) {
            ret = true;
            break;
        }
    }
    return ret;
}


bool Node::addChild (Node& c) {
    if (hasChild(c) == false) {
        children.push_back(&c);
        c.setParent(*this);
        return true;
    } else {
        return false;
    }
}


bool Node::removeChild (Node& c) {
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


Node * Node::getChild (int c) {
    return children.at(c);
}


std::string Node::getName () {
    return name;
}


void Node::setName (std::string s) {
    name = s;
}


std::string Node::getComment () {
    return comment;
}

void Node::setComment (std::string s) {
    comment = s;
}


// since nexus writer uses this, using stricter nexus punctuation
std::string Node::getNewick (bool bl) {
    std::string ret = "";
    for (int i=0; i < this->getChildCount(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret = ret+this->getChild(i)->getNewick(bl);
        if (bl == true) {
            //std::ostringstream o;
            ////20 is what you get from raxml
            //o.setf(ios::fixed, std::ios::floatfield);
            //o << std::setprecision(20) << this->getChild(i)->getBL();
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
        // newick punct is a subset of Nexus, so labels will be safe
        std::string compliantName = get_valid_nexus_label(name);
        ret += compliantName;
    }
    return ret;
}


// atm returns both 1) knuckles and 2) root edges
// neither of these is likely wanted
// root edge could be solved by passing in the correct node initially, but MRCA is expensive
// recursive
std::string Node::getPaintedNewick (bool bl) {
    std::string ret = "";
    std::vector<int> paintedchildren;
    for (int i=0; i < this->getChildCount(); i++) {
        if (this->getChild(i)->getPainted() == true) {
            paintedchildren.push_back(i);
        }
    }
    /*
    string ndname = this->getName();
    if (ndname != "") {
            std::cout << "Dealing with node '" << ndname << "'." << std::endl;
        }
    std::cout << "   dealing with " << paintedchildren.size() << " children." << std::endl;
    
    if (paintedchildren.size() == 1) {
        std::cout << "looks like '" << ndname << "' here might be a knuckle!" << std::endl;
        bool done = false;
        double el = this->getChild(paintedchildren[0])->getBL();
        while (!done) {
            Node * cnd = this->getChild(paintedchildren[0]);
        }
    }
    */
    for (unsigned int i=0; i < paintedchildren.size(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret += this->getChild(paintedchildren[i])->getPaintedNewick(bl);
        if (bl == true) {
            std::ostringstream o;
            // 20 is what you get from raxml
            o << std::setprecision(20) << this->getChild(paintedchildren[i])->getBL();
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
std::string Node::getNewick (bool bl, std::string obj) {
    std::string ret = "";
    for (int i=0; i < this->getChildCount(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret = ret+this->getChild(i)->getNewick(bl, obj);
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


Node * Node::getParent () {
    return parent;
}


int Node::getChildCount () {
    return children.size();
}


void Node::assocObject (std::string name, NodeObject& obj) {
    if (assoc.count(name) > 0) {
        delete assoc[name];
    }
    assoc[name] = obj.clone();
}


void Node::assocDoubleVector (std::string name, std::vector<Superdouble>& obj) {
    if (assocDV.count(name) > 0) {
        assocDV.erase(name);
    }
    std::vector<Superdouble> tvec (obj.size());
    for (unsigned int i=0; i < obj.size(); i++) {
        tvec[i] = obj[i];
    }
    assocDV[name] = tvec;
}


std::vector<Superdouble> * Node::getDoubleVector (std::string name) {
    return &assocDV[name];
}


void Node::deleteDoubleVector (std::string name) {
    if (assocDV.count(name) > 0) {
        assocDV.erase(name);
    }
}


/*
 * gets the number of leaves from this node
 */
int Node::get_num_leaves () {
    int retnum = 0;
    std::stack<Node*> nodes;
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
std::vector<Node*> Node::get_leaves () {
    std::vector<Node*> retnodes;
    std::stack<Node*> nodes;
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
std::set<Node*> Node::get_leaves_set () {
    std::set<Node*> retnodes;
    std::stack<Node*> nodes;
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


std::set<std::string> Node::get_leave_names_set () {
    std::stack<Node*> nodes;
    std::set<std::string> names;
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


std::vector<std::string> Node::get_leave_names () {
    std::stack<Node*> nodes;
    std::vector<std::string> names;
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
 * tree.getRoot()->assocObject("test", sno);
 * std::cout << *((StringNodeObject*) (tree.getRoot()->getObject("test"))) << std::endl;
 *
 * and the vector like
 * VectorNodeObject<int> vno;
 * vno.push_back(1);vno.push_back(2);
 * tree.getRoot()->assocObject("testvno", vno);
 * std::cout << ((VectorNodeObject<int> *) (tree.getRoot()->getObject("testvno")))->at(0) << std::endl;
 */
NodeObject  * Node::getObject (std::string name) {
    return assoc[name];
}


void Node::initSegVector () {
    segs = new std::vector<BranchSegment> ();
}


std::vector<BranchSegment> * Node::getSegVector () {
    return segs;
}


void Node::deleteSegVector () {
    delete segs;
}


/*
 * delete the node
 */
Node::~Node () {
    std::map<std::string, NodeObject *>::iterator it;
    for (it = assoc.begin() ; it != assoc.end(); it++) {
        delete assoc[it->first];
    }
}
