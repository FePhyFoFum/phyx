#ifndef PX_NODE_H
#define PX_NODE_H

#include <map>
#include <set>
#include <string>

#include "branch_segment.h"
#include "vector_node_object.h"
#include "superdouble.h"

class Node {
private:
    double BL_; // branch length, points to parent
    double height_; // could be from tip or from root. need to be consistent
    double depth_; // not being used. should be
    unsigned int number_;
    std::string name_;
    Node * parent_;
    // children are _immediate_ children, not all descendants
    std::vector<Node *> children_;
    std::map<std::string, NodeObject *> assoc_;
    std::map<std::string, std::vector<Superdouble> > assocDV_;
    std::vector<BranchSegment> * segs_;
    std::string comment_;
    bool painted_;

public:
    Node ();
    explicit Node (Node * inparent);
    // not used
    Node (double bl, unsigned int innumber, std::string inname, Node * inparent);
    
    int get_num_leaves ();
    std::vector<Node*> get_leaves ();
    std::vector<std::string> get_leave_names ();
    std::set<std::string> get_leave_names_set ();
    std::set<Node *> get_leaves_set ();
    
    std::vector<Node*> getChildren () const;
    bool isExternal ();
    bool isInternal ();
    bool isRoot ();
    bool isKnuckle ();
    bool hasParent ();
    void setParent (Node& p);
    unsigned int getNumber () const;
    void setNumber (unsigned int n);
    bool getPainted () const;
    void setPainted (bool p);
    double getBL () const;
    void setBL (double bl);
    double getHeight () const;
    void setHeight (double he);
    double getDepth () const;
    void setDepth (double de);
    bool hasChild (Node& test);
    bool addChild (Node& c);
    bool removeChild (Node& c);
    Node * getChild (unsigned int c);
    std::string getName () const;
    std::string getComment () const;
    void setName (std::string s);
    void setComment (std::string s);
    std::string getNewick (bool bl);
    std::string getNewick (bool bl, const std::string& obj);
    std::string getPaintedNewick (bool bl);
    Node * getParent () const;
    unsigned int getChildCount ();
    void assocObject (const std::string& name, NodeObject& obj);
    void assocDoubleVector (const std::string& name, std::vector<Superdouble>& obj);
    std::vector<Superdouble> * getDoubleVector (const std::string& name);
    void deleteDoubleVector (const std::string& name);
    NodeObject * getObject (const std::string& name);
    void initSegVector ();
    std::vector<BranchSegment> * getSegVector () const;
    void deleteSegVector ();
    
    VectorNodeObject<Superdouble> seg_sp_stoch_map_revB_time; //segment specific rev B, combining the tempA and the ENLT
    VectorNodeObject<Superdouble> seg_sp_stoch_map_revB_number; //segment specific rev B, combining the tempA and the ENLT
    ~Node ();
};

#endif /* PX_NODE_H */
