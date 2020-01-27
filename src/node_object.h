#ifndef _NODE_OBJECT_H_
#define _NODE_OBJECT_H_


class NodeObject {

public:
    NodeObject() {}

    virtual ~NodeObject() {}

public:
    virtual NodeObject * clone() const = 0;
};

#endif /* _NODE_OBJECT_H_ */
