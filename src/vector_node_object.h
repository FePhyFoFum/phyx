/*
 * node_object.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef VECTOR_NODE_OBJECT_H_
#define VECTOR_NODE_OBJECT_H_

#include "node_object.h"
#include <vector>
using namespace std;

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

#endif /* NODE_OBJECT_H_ */
