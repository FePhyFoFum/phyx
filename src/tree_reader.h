#ifndef _TREE_READER_H_
#define _TREE_READER_H_

#include <string>
#include <map>
#include <iostream>

#include "tree.h"


class TreeReader {
public:
    TreeReader ();
    Tree * readTree (const std::string& trees);
};


// hrm why are these not member functions?
Tree * read_tree_string (std::string pb);
int test_tree_filetype (const std::string&);
int test_tree_filetype_stream (std::istream& stri, std::string& retstring);
bool get_nexus_translation_table (std::istream& stri,
        std::map<std::string, std::string> * trans, std::string * retstring);
Tree * read_next_tree_from_stream_nexus (std::istream& stri, std::string& retstring,
        bool ttexists, std::map<std::string, std::string> * trans, bool * going);
Tree * read_next_tree_from_stream_newick (std::istream& stri, std::string& retstring,
        bool * going);

#endif /* _TREE_READER_H_ */
