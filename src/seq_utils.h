#ifndef _SEQ_UTILS_H_
#define _SEQ_UTILS_H_

#include <map>
#include <set>

using namespace std;

#include "sequence.h"

string guess_alignment_type (string & sequence);
char get_dna_from_pos(set<int> ins);
set<int> get_dna_pos(char);
string consensus_seq(vector<Sequence> &, string & alpha);
char single_dna_complement(char inc);
void write_phylip_alignment(vector<Sequence> & seqs, bool const& uppercase, ostream * ostr);
void write_nexus_alignment(vector<Sequence> & seqs, bool const& uppercase, ostream * ostr);
void populate_codon_list(vector<string> * codon_list);
void populate_map_codon_dict(map<string, string> * codon_dict);
void populate_map_codon_indices(map<string, vector<int> > * codon_position);
void create_vector_seq_codon_state_reconstructor(vector<Sequence> & origseqs,
    vector<Sequence> & sr_seqs, int site, map<string,vector<int> > & codon_pos);

#endif /* _SEQ_UTILS_H_ */
