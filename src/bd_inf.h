#ifndef _BD_INF_H_
#define _BD_INF_H_

#include <string>

using namespace std;

#include "tree.h"

class BirthDeathInference {
private:
    string model;
    Tree* tree;
    double treelength;
    double nintnodes;
    double nspeciation;
    double lambda;
    double epsilon;
    
    



public:
    BirthDeathInference(string const& model);
    
};

#endif /* _BD_INF_H_ */
