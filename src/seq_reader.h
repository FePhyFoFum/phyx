#ifndef SEQ_READER_H_
#define SEQ_READER_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "sequence.h"

using namespace std;

int test_seq_filetype(string);
int test_seq_filetype_stream(istream &, string &);
int test_char_filetype_stream(istream & stri,string & retstring);
bool read_next_seq_from_stream(istream & stri, int ftype, string & retstring, Sequence & seq);
bool read_next_seq_char_from_stream(istream & stri, int ftype, string & retstring, Sequence & seq);
bool read_fasta_file(string,vector<Sequence>&);
bool read_phylip_file(string,vector<Sequence>&);
bool read_phylip_file_strec(string,vector<Sequence>&);

#endif
