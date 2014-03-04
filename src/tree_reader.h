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
#include <map>

using namespace std;

#include "node.h"
#include "tree.h"

class TreeReader{
public:
	TreeReader();
	Tree * readTree(string trees);
};

Tree * read_tree_string(string trees);
int test_tree_filetype(string filen);
int test_tree_filetype_stream(istream & stri, string & retstring);
bool get_nexus_translation_table(istream & stri, map<int,string> * trans);
bool read_next_tree_from_stream_nexus(istream & stri, string & retstring, Tree & tree);
bool read_next_tree_from_stream_newick(istream & stri, string & retstring, Tree & tree);

#endif /* TREE_READER_H_ */
