/*
 * seqgen.h
 *
 *  Created on: Jun 23, 2015
 *      Author: joe
 */

#ifndef SEQGEN_H_
#define SEQGEN_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <map>
#include <iterator>

#include "seqgen.h"
#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "node.h"
#include "tree.h"
#include "tree_reader.h"


class SEQGEN {

private:

    vector<Node *> nodes;


public:
	SEQGEN();
	void TakeInTree(vector< vector<double> >&, Tree *, int, vector<double>&);
	Node * PreOrder(int);
	virtual ~SEQGEN();
};

#endif /* SEQGEN_H_ */
