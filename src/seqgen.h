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
#include <random>

//#include "seqgen.h"
//#include "sequence.h"
//#include "seq_reader.h"
//#include "utils.h"
//#include "node.h"
#include "tree.h"
//#include "tree_reader.h"


class SequenceGenerator {

private:
    // constant values
    int seqlen;
    vector <double> basefreqs;
    vector < vector <double> > rmatrix;
    Tree * tree;
    int nreps;
    int seed;
    string Ancestor;
    static map<char, int> nucMap;
    bool showancs;
    
    map <Node *, string> seqs;
    
    void EvoSim (Tree * tree);
    vector < vector <double> > CalQ ();
    vector < vector <double> > PCalq (vector < vector <double> > QMatrix, float br);
    string SeqSim (string const& anc, vector < vector <double> >& Matrix);
    void randDNA ();
    //void TakeInTree ();
    mt19937 generator;
    float get_uniform_random_deviate ();
//    vector<Node *> nodes; // unused


public:
    SequenceGenerator ();
    SequenceGenerator (int const &seqlenth, vector <double> const& basefreq,
        vector < vector<double> >& rmatrix, Tree * tree, bool const& showancs, 
        int const& nreps, int const & seed);
    //Node * PreOrder(int); // not used
    
    
    //virtual ~SEQGEN();
};

#endif /* _SEQGEN_H_ */
