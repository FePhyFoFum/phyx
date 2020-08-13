#ifndef _SEQ_READER_H_
#define _SEQ_READER_H_

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
void get_nexus_dimensions (std::istream& stri, int& num_taxa, int& numChar, bool& interleave);
void get_nexus_alignment_properties (std::istream& stri, int& num_taxa, int& numChar,
        bool& interleave, std::string& alpha_name, std::string& symbols, char& gap, char& missing);
std::vector<Sequence> read_interleaved_nexus (std::istream& stri, int num_taxa, int num_char);
void get_phylip_dimensions (std::string head, int& num_taxa, int& numChar);
bool is_complicated_phylip (std::istream& pios, const int& num_char);
void get_phylip_format (std::istream& pios, const unsigned int& num_taxa, const unsigned int& numChar,
        bool& interleaved, bool& spaces, bool& multiline);
std::vector<Sequence> read_phylip (std::istream& pios, const int& num_taxa, const int& numChar);
std::vector<Sequence> ingest_alignment (std::istream* pios, std::string& alphaName);

// deprecated
//int test_seq_filetype (std::string filen);
//bool read_fasta_file (std::string filen, std::vector<Sequence>& seqs);
//bool read_phylip_file (std::string filen, std::vector<Sequence>& seqs);
//bool read_phylip_file_strec (std::string filen, std::vector<Sequence>& seqs);//bool read_nexus_seqs_file (std::string filen, std::vector<Sequence>& seqs);
//void get_nexus_dimensions_file (std::string& filen, int& num_taxa, int& numChar, bool& interleave);
//std::vector<Sequence> read_interleaved_nexus_file (std::string filen, int num_taxa, int num_char);

#endif /* _SEQ_READER_H_ */
