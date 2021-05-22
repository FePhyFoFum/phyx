#ifndef PX__SEQ_READER_H
#define PX__SEQ_READER_H

#include <string>
#include <vector>
#include <iostream>

class Sequence; // forward declaration

std::string get_filetype_string (const int& ft);
int test_seq_filetype_stream (std::istream& stri, std::string& retstring);
int test_char_filetype_stream (std::istream& stri, std::string& retstring);
bool read_next_seq_from_stream (std::istream& stri, int ftype, std::string& retstring,
    Sequence& seq);
bool read_next_seq_char_from_stream (std::istream& stri, int ftype,
    std::string& retstring, Sequence& seq);
void get_nexus_dimensions (std::istream& stri, int& num_taxa, int& num_char, bool& interleave);
void get_nexus_alignment_properties (std::istream& stri, int& num_taxa, int& num_char,
        bool& interleave, std::string& alpha_name, std::string& symbols, char& gap, char& missing);
std::vector<Sequence> read_interleaved_nexus (std::istream& stri, int num_taxa, int num_char);
void get_phylip_dimensions (const std::string& head, int& num_taxa, int& num_char);
bool is_complicated_phylip (std::istream& pios, const int& num_char);
void get_phylip_format (std::istream& pios, const unsigned int& num_char,
        bool& interleaved, bool& spaces, bool& multiline);
std::vector<Sequence> read_phylip (std::istream& pios, const int& num_taxa, const int& num_char);
std::vector<Sequence> ingest_alignment (std::istream* pios, std::string& alphaName);

#endif /* PX__SEQ_READER_H */
