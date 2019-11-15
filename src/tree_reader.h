/*
 * tree_reader.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef _TREE_READER_H_
#define _TREE_READER_H_

#include <string>
#include <map>
#include <iostream>

#include "node.h"
#include "tree.h"

class TreeReader {
public:
    TreeReader ();
    Tree * readTree (std::string trees);
};

Tree * read_tree_string (std::string trees);
int test_tree_filetype (std::string filen);
int test_tree_filetype_stream (std::istream& stri, std::string& retstring);
bool get_nexus_translation_table (std::istream& stri,
        std::map<std::string, std::string> * trans, std::string * retstring);
bool check_nexus_comment (std::string line);
void process_nexus_comment (std::istream& stri, std::string& tline);
Tree * read_next_tree_from_stream_nexus (std::istream& stri, std::string& retstring,
        bool ttexists, std::map<std::string, std::string> * trans, bool * going);
Tree * read_next_tree_from_stream_newick (std::istream& stri, std::string& retstring,
        bool * going);

#endif /* _TREE_READER_H_ */
