/*
 * node_object.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef PX_STRING_NODE_OBJECT_H
#define PX_STRING_NODE_OBJECT_H

#include <string>

#include "node_object.h"

class StringNodeObject: public std::string, public NodeObject {
public:
    StringNodeObject(const char * value): std::string(value) {}
    StringNodeObject(const std::string& value): std::string(value) {}
    virtual ~StringNodeObject() {}

public:
    StringNodeObject * clone() const { return new StringNodeObject(*this); }
};

#endif /* PX_STRING_NODE_OBJECT_H */
