#ifndef PX_VECTOR_NODE_OBJECT_H
#define PX_VECTOR_NODE_OBJECT_H

#include <vector>

#include "node_object.h"


template<class T>
class VectorNodeObject: public std::vector<T>, public NodeObject {
public:
    VectorNodeObject():
        std::vector<T> () {}
    
    VectorNodeObject(typename std::vector<T>::size_type num, const T& val = T() ):
        std::vector<T> (num, val) {}
    
    VectorNodeObject(typename std::vector<T>::iterator start, typename std::vector<T>::iterator end):
        std::vector<T> (start, end) {}
    
    virtual ~VectorNodeObject () {}
    
public:
    VectorNodeObject<T> * clone () const { return new VectorNodeObject<T>(*this); }
};

#endif /* PX_VECTOR_NODE_OBJECT_H */
