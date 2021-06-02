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


Node::Node ():BL_(0.0), height_(0.0), depth_(0.0), number_(0), name_(""), parent_(nullptr),
        children_(std::vector<Node *> ()), assoc_(std::map<std::string, NodeObject *>()),
        assocDV_(std::map<std::string, std::vector<Superdouble> >()), segs_(nullptr),
        comment_(""), painted_(false) {

}


Node::Node (Node * inparent):BL_(0.0), height_(0.0), depth_(0.0), number_(0), name_(""),
        parent_(inparent), children_(std::vector<Node *> ()), assoc_(std::map<std::string,
        NodeObject *>()), assocDV_(std::map<std::string, std::vector<Superdouble> >()),
        segs_(nullptr), comment_(""), painted_(false) {

}


Node::Node (double bl, int innumber, std::string inname, Node * inparent):BL_(bl),
        height_(0.0), depth_(0.0), number_(innumber), name_(std::move(inname)),
        parent_(inparent), children_(std::vector<Node *> ()),
        assoc_(std::map<std::string, NodeObject *>()),
        assocDV_(std::map<std::string, std::vector<Superdouble> >()), segs_(nullptr),
        comment_(""), painted_(false) {

}


std::vector<Node*> Node::getChildren () {
    return children_;
}


bool Node::isExternal () {
    bool ret = false;
    if (children_.empty()) {
        ret = true;
    }
    return ret;
}


bool Node::isInternal () {
    bool ret = false;
    if (!children_.empty()) {
        ret = true;
    }
    return ret;
}


bool Node::isRoot () {
    bool ret = false;
    if (parent_ == nullptr) {
        ret = true;
    }
    return ret;
}


// is the internal node of two degree (one parent and one child?)
// neither tips nor the root can be a knuckle
bool Node::isKnuckle () {
    bool ret = false;
    if (!isRoot() && !isExternal()) {
        if (children_.size() == 1) {
            ret = true;
        }
    }
    return ret;
}


bool Node::hasParent () {
    bool ret = true;
    if (parent_ == nullptr) {
        ret = false;
    }
    return ret;
}


void Node::setParent (Node& p) {
    parent_ = &p;
}


int Node::getNumber () {
    return number_;
}


void Node::setNumber (int n) {
    number_ = n;
}


bool Node::getPainted () {
    return painted_;
}


void Node::setPainted (bool p) {
    painted_ = p;
}


double Node::getBL () {
    return BL_;
}


void Node::setBL (double bl) {
    BL_ = bl;
}


double Node::getHeight () {
    return height_;
}


void Node::setHeight (double he) {
    height_ = he;
}


double Node::getDepth () {
    return depth_;
}


void Node::setDepth (double de) {
    depth_ = de;
}


bool Node::hasChild (Node& test) {
    bool ret = false;
    for (auto & chi : children_) {
        if (chi == &test) {
            ret = true;
            break;
        }
    }
    return ret;
}


bool Node::addChild (Node& c) {
    bool ret = false;
    if (!hasChild(c)) {
        children_.push_back(&c);
        c.setParent(*this);
        ret = true;
    }
    return ret;
}


bool Node::removeChild (Node& c) {
    bool ret = false;
    if (hasChild(c)) {
        for (unsigned int i = 0; i < children_.size(); i++) {
            if (children_.at(i) == &c) {
                children_.erase(children_.begin()+i);
                break;
            }
        }
        ret = true;
    }
    return ret;
}


Node * Node::getChild (int c) {
    return children_.at(static_cast<unsigned long>(c));
}


std::string Node::getName () {
    return name_;
}


void Node::setName (std::string s) {
    name_ = std::move(s);
}


std::string Node::getComment () {
    return comment_;
}

void Node::setComment (std::string s) {
    comment_ = std::move(s);
}


// since nexus writer uses this, using stricter nexus punctuation
std::string Node::getNewick (bool bl) {
    std::string ret;
    for (int i = 0; i < this->getChildCount(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret += this->getChild(i)->getNewick(bl);
        if (bl) {
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
    if (!name_.empty()) {
        // newick punct is a subset of Nexus, so labels will be safe
        std::string compliantName = get_valid_nexus_label(name_);
        ret += compliantName;
    }
    return ret;
}


// atm returns both 1) knuckles and 2) root edges
// neither of these is likely wanted
// root edge could be solved by passing in the correct node initially, but MRCA is expensive
// recursive
std::string Node::getPaintedNewick (bool bl) {
    std::string ret;
    std::vector<int> paintedchildren;
    for (int i = 0; i < this->getChildCount(); i++) {
        if (this->getChild(i)->getPainted()) {
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
    for (unsigned int i = 0; i < paintedchildren.size(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret += this->getChild(paintedchildren[i])->getPaintedNewick(bl);
        if (bl) {
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
    if (!this->getName().empty()) {
        ret += this->getName();
    }
    return ret;
}


/*
 * should be returning the stringnodeobjects as the names for internal
 * nodes with the [&obj=string]
 * needs to be a string in setObject
 */
std::string Node::getNewick (bool bl, const std::string& obj) {
    std::string ret;
    for (int i = 0; i < this->getChildCount(); i++) {
        if (i == 0) {
            ret += "(";
        }
        ret += this->getChild(i)->getNewick(bl, obj);
        if (bl) {
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
    if (!this->name_.empty()) {
        ret = ret + this->name_;
    }
    if (this->getObject(obj) != nullptr) {
        std::ostringstream o;
        o << (*(dynamic_cast<StringNodeObject*>(this->getObject(obj))));
        ret += "[&"+obj+"="+o.str()+"]";
    }
    return ret;
}


Node * Node::getParent () {
    return parent_;
}


int Node::getChildCount () {
    return static_cast<int>(children_.size());
}


void Node::assocObject (const std::string& name, NodeObject& obj) {
    if (assoc_.count(name) > 0) {
        delete assoc_[name];
    }
    assoc_[name] = obj.clone();
}


void Node::assocDoubleVector (const std::string& name, std::vector<Superdouble>& obj) {
    if (assocDV_.count(name) > 0) {
        assocDV_.erase(name);
    }
    std::vector<Superdouble> tvec (obj.size());
    for (unsigned int i = 0; i < obj.size(); i++) {
        tvec[i] = obj[i];
    }
    assocDV_[name] = tvec;
}


std::vector<Superdouble> * Node::getDoubleVector (const std::string& name) {
    return &assocDV_[name];
}


void Node::deleteDoubleVector (const std::string& name) {
    if (assocDV_.count(name) > 0) {
        assocDV_.erase(name);
    }
}


/*
 * gets the number of leaves from this node
 */
int Node::get_num_leaves () {
    int retnum = 0;
    std::stack<Node*> nodes;
    nodes.push(this);
    while (!nodes.empty()) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i = 0; i < nd->getChildCount(); i++) {
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
    while (!nodes.empty()) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i = 0; i < nd->getChildCount(); i++) {
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
    while (!nodes.empty()) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i = 0; i < nd->getChildCount(); i++) {
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
    while (!nodes.empty()) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i = 0; i < nd->getChildCount(); i++) {
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
    while (!nodes.empty()) {
        Node * nd = nodes.top();
        nodes.pop();
        if (nd->getChildCount() > 0) {
            for (int i = 0; i < nd->getChildCount(); i++) {
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
NodeObject  * Node::getObject (const std::string& name) {
    return assoc_[name];
}


void Node::initSegVector () {
    segs_ = new std::vector<BranchSegment> ();
}


std::vector<BranchSegment> * Node::getSegVector () {
    return segs_;
}


void Node::deleteSegVector () {
    delete segs_;
}


/*
 * delete the node
 */
Node::~Node () {
    std::map<std::string, NodeObject *>::iterator it;
    for (it = assoc_.begin(); it != assoc_.end(); ++it) {
        delete assoc_[it->first];
    }
}
