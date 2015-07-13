/*
 * seqgen.h
 *
 *  Created on: Jun 23, 2015
 *      Author: joe
 */

#ifndef _SEQGEN_H_
#define _SEQGEN_H_

//#include <iostream>
//#include <string>
//#include <fstream>
//#include <vector>
//#include <sstream>
//#include <iterator>
//#include <algorithm>
//#include <map>
//#include <iterator>

//#include "seqgen.h"
//#include "sequence.h"
//#include "seq_reader.h"
//#include "utils.h"
//#include "node.h"
#include "tree.h"
//#include "tree_reader.h"


class SEQGEN {

private:
    // constant values
    int seqlen;
    vector <double> basefreqs;
    vector < vector <double> > rmatrix;
    Tree * tree;
    int nreps;
    
    void EvoSim(Tree * tree, string Ancestor);
    vector < vector <double> > CalQ();
    vector < vector <double> > PCalq(vector < vector <double> > QMatrix, float br);
    void SeqSim(string& Ancestor, vector < vector <double> >& Matrix);
    void randDNA (string& Ancestor);
    void TakeInTree();
//    vector<Node *> nodes; // unused


public:
    SEQGEN();
    SEQGEN(int const &seqlenth, vector <double> const& basefreq,
        vector < vector<double> >& rmatrix, Tree * tree, int const& nreps);
    //Node * PreOrder(int);
    
    
    //virtual ~SEQGEN();
};

#endif /* _SEQGEN_H_ */
