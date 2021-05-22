#ifndef PX__NODE_OBJECT_H
#define PX__NODE_OBJECT_H


class NodeObject {

public:
    NodeObject() {}

    virtual ~NodeObject() {}

public:
    virtual NodeObject * clone() const = 0;
};

#endif /* PX__NODE_OBJECT_H */
