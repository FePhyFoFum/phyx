#ifndef _NODE_H_
#define _NODE_H_

#include <map>
#include <set>
#include <string>

#include "branch_segment.h"
#include "vector_node_object.h"
#include "superdouble.h"

class Node {
private:
    double BL; // branch length, points to parent
    double height; // could be from tip or from root
    double depth; // not being used
    int number;
    std::string name;
    Node * parent;
    std::vector<Node *> children;
    std::map<std::string, NodeObject *> assoc;
    std::map<std::string, std::vector<Superdouble> > assocDV;
    std::vector<BranchSegment> * segs;
    std::string comment;
    bool painted;

public:
    Node ();
    Node (Node * parent);
    Node (double bl, int number, std::string name, Node * parent);
    
    int get_num_leaves ();
    std::vector<Node*> get_leaves ();
    std::vector<std::string> get_leave_names ();
    std::set<std::string> get_leave_names_set ();
    std::set<Node *> get_leaves_set ();
    
    std::vector<Node*> getChildren ();
    bool isExternal ();
    bool isInternal ();
    bool isRoot ();
    bool isKnuckle ();
    bool hasParent ();
    void setParent (Node& p);
    int getNumber ();
    void setNumber (int n);
    bool getPainted ();
    void setPainted (bool p);
    double getBL ();
    void setBL (double bl);
    double getHeight ();
    void setHeight (double he);
    double getDepth ();
    void setDepth (double de);
    bool hasChild (Node& test);
    bool addChild (Node& c);
    bool removeChild (Node& c);
    Node * getChild (int c);
    std::string getName ();
    std::string getComment ();
    void setName (std::string s);
    void setComment (std::string s);
    std::string getNewick (bool bl);
    std::string getNewick (bool bl, std::string obj);
    std::string getPaintedNewick (bool bl);
    Node * getParent ();
    int getChildCount ();
    void assocObject (std::string name, NodeObject& obj);
    void assocDoubleVector (std::string name, std::vector<Superdouble>& obj);
    std::vector<Superdouble> * getDoubleVector (std::string name);
    void deleteDoubleVector (std::string name);
    NodeObject * getObject (std::string name);
    void initSegVector ();
    std::vector<BranchSegment> * getSegVector ();
    void deleteSegVector ();
    
    VectorNodeObject<Superdouble> seg_sp_stoch_map_revB_time; //segment specific rev B, combining the tempA and the ENLT
    VectorNodeObject<Superdouble> seg_sp_stoch_map_revB_number; //segment specific rev B, combining the tempA and the ENLT
    ~Node ();
};

#endif /* _NODE_H_ */
