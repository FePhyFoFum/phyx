/*
 * node_object.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef _NODE_OBJECT_H_
#define _NODE_OBJECT_H_

using namespace std;

class NodeObject{

public:
    NodeObject() {}

    virtual ~NodeObject() {}

public:
    virtual NodeObject * clone() const = 0;
};

#endif /* _NODE_OBJECT_H_ */
