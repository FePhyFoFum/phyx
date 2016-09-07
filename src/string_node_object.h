/*
 * node_object.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef _STRING_NODE_OBJECT_H_
#define _STRING_NODE_OBJECT_H_

using namespace std;

#include "node_object.h"

class StringNodeObject: public string, public NodeObject{
public:
    StringNodeObject(const char * value): string(value) {}
    StringNodeObject(const string & value): string(value) {}
    virtual ~StringNodeObject() {}

public:
    StringNodeObject * clone() const { return new StringNodeObject(*this); }
};

#endif /* _STRING_NODE_OBJECT_H_ */
