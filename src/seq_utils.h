#ifndef PX_SEQ_UTILS_H
#define PX_SEQ_UTILS_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>

class Sequence; // forward declarations


std::string guess_alignment_type (std::string& sequence);
char get_dna_from_pos (const std::set<int>& ins);
char get_prot_char (const std::set<char>& inc);
std::set<int> get_dna_pos (char);
std::string consensus_seq (std::vector<Sequence>&, std::string& alpha);
char single_dna_complement (char inc);
void write_phylip_alignment (std::vector<Sequence>& seqs, const bool& uppercase, std::ostream * ostr);
void write_nexus_alignment (std::vector<Sequence>& seqs, const bool& uppercase, std::ostream * ostr);
std::vector<std::string> collect_names (const std::vector<Sequence>& algnmnt);
void populate_codon_list (std::vector<std::string> * codon_list);
void populate_map_codon_dict (std::map<std::string, std::string> * codon_dict);
void populate_map_codon_indices (std::map<std::string, std::vector<int> > * codon_position);
void create_vector_seq_codon_state_reconstructor (std::vector<Sequence>& origseqs,
    std::vector<Sequence>& sr_seqs, int site, std::map<std::string, std::vector<int> >& codon_pos);
bool check_binary_sequence (const std::string& seq);
std::string get_alphabet_from_sequence (const std::string& instr);
bool is_dna_char (char& residue);
bool is_prot_char (char& residue);
int count_dna_chars (const std::string& str);
bool is_aligned (const std::vector<Sequence>& seqs);
bool is_codon_alignment (const std::vector<Sequence>& seqs);

#endif /* PX_SEQ_UTILS_H */
