#ifndef PX_TREE_READER_H
#define PX_TREE_READER_H

#include <string>
#include <map>
#include <iostream>

#include "tree.h"


class TreeReader {
public:
    TreeReader ();
    Tree * readTree (const std::string& pb);
};


// hrm why are these not member functions?
int test_tree_filetype (const std::string&);
int test_tree_filetype_stream (std::istream& stri, std::string& retstring);
bool get_nexus_translation_table (std::istream& stri,
        std::map<std::string, std::string> * trans, std::string * retstring);
Tree * read_next_tree_from_stream_nexus (std::istream& stri, std::string& retstring,
        bool ttexists, std::map<std::string, std::string> * trans, bool * going);
Tree * read_next_tree_from_stream_newick (std::istream& stri, std::string& retstring,
        bool * going);

#endif /* PX_TREE_READER_H */
