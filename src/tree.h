#ifndef PX_TREE_H
#define PX_TREE_H

#include <string>
#include <vector>
#include <map>

#include "node.h"


class Tree {
private:
    Node * root_;
    std::vector<Node *> nodes_;
    std::vector<Node *> internal_nodes_;
    std::vector<Node *> external_nodes_;
    std::map<std::string, Node*> name_node_map_;
    unsigned int internal_node_count_;
    unsigned int external_node_count_;
    bool edge_lengths_;
    bool node_annotations_;
    bool internal_node_names_;
    
    void processReRoot (Node * node);
    void exchangeInfo (Node * node1, Node * node2);
    void postOrderProcessRoot (Node * node);
    Node * getMRCATraverse (Node * curn1, Node * curn2);
    void setHeightFromRootToNode (Node& inNode, double newHeight);
    double getGreatestDistance (Node * inNode);
    
public:
    Tree ();
    explicit Tree (Node * inroot);
    Tree * clone () { return new Tree(*this); }
    
    void addExternalNode (Node * tn);
    void addInternalNode (Node * tn);
    void pruneExternalNode (Node * node);
    void pruneInternalNode (Node * node);
    Node * getExternalNode (unsigned int num);
    Node * getExternalNode (const std::string& name);
    Node * getInternalNode (unsigned int num);
    Node * getInternalNode (std::string& name);
    Node * getNode (unsigned int num);
    Node * getNode (std::string& name);
    unsigned int getNodeCount () const;
    unsigned int getExtantNodeCount ();
    unsigned int getExternalNodeCount () const;
    unsigned int getInternalNodeCount () const;
    Node * getRoot () const;
    void setRoot (Node * inroot);
    void setEdgeLengthsPresent (bool res);
    bool hasEdgeLengths () const;
    void setNodeAnnotationsPresent (bool res);
    bool hasNodeAnnotations () const;
    void setNodeNamesPresent (bool res);
    bool hasNodeNames () const;
    void unRoot ();
    bool reRoot (Node * inroot);
    void duplicateRootSupport ();
    void tritomyRoot (Node * toberoot);
    Node * getMRCA (std::vector<std::string> innodes);
    Node * getMRCA (std::vector<Node *> innodes);
    Node * getInternalMRCA (std::vector<std::string>& innodes);
    void processRoot ();
    void exchangeNodes (Node * node1, Node * node2);
    void removeRootEdge ();
    void setHeightFromRootToNodes ();
    void setHeightFromTipToNodes ();
    
    ~Tree ();
};

#endif /* PX_TREE_H */
