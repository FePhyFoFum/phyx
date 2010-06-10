/*
 * tree_reader.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef TREE_READER_H_
#define TREE_READER_H_

#include <string>
#include <vector>

using namespace std;

#include "node.h"
#include "tree.h"

class TreeReader{
public:
	TreeReader();
	Tree * readTree(string trees);
};

#endif /* TREE_READER_H_ */
