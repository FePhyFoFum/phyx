/*
 * node_object.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef _VECTOR_NODE_OBJECT_H_
#define _VECTOR_NODE_OBJECT_H_

#include <vector>

using namespace std;

#include "node_object.h"

template<class T>
class VectorNodeObject: public vector<T>, public NodeObject{
public:
    VectorNodeObject():
        vector<T>() {}

    VectorNodeObject(typename vector<T>::size_type num, const T& val = T() ):
        vector<T>(num, val) {}

    VectorNodeObject(typename vector<T>::iterator start, typename vector<T>::iterator end):
        vector<T>(start, end) {}

    virtual ~VectorNodeObject() {}

public:

    VectorNodeObject<T> * clone() const { return new VectorNodeObject<T>(*this); }
};

#endif /* _VECTOR_NODE_OBJECT_H_ */
