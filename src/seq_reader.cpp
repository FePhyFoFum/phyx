#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <sstream>

#include "sequence.h"
#include "seq_reader.h"
#include "seq_utils.h"
#include "utils.h"


// for printing purposes (and to remind programmers what ints refer to)
std::map<int, std::string> filetype_map = {
   {0, "nexus"},
   {1, "phylip"},
   {2, "fasta"},
   {3, "fastq"}
};


std::string get_filetype_string (const int& ft) {
    std::string ftype = filetype_map[ft];
    return ftype;
}

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta), @ (fastq)
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
*/
int test_seq_filetype_stream(std::istream& stri, std::string& retstring) {
    if (!getline_safe(stri, retstring)) {
        std::cout << "ERROR: end of file too soon" << std::endl;
    }
    int ret = 666; // if you get 666, there is no filetype set
    if (retstring[0] == '#') {
        ret = 0;
    } else if (retstring[0] == '>') {
        ret = 2;
    } else if (retstring[0] == '@') {
        ret = 3;
    } else {
        std::vector<std::string> tokens;
        std::string del(" \t");
        tokenize(retstring, tokens, del);
        if (tokens.size() > 1) {
            trim_spaces(tokens[0]);
            if (is_number(tokens[0])) {
                ret = 1;
            }
        }
    }
    return ret;
}


/*
 * returns the next string in the getline if there is one
 * TODO: nexus interleaved is not going to work here (read_interleaved_nexus does this)
 * TODO: - skip Nexus comment lines - done
 *       - check if label is quoted. if so, need to skip any whitespace that might be present
*/
bool read_next_seq_from_stream (std::istream & stri, int ftype, std::string& retstring, Sequence& seq) {
    std::string tline;
    if (ftype == 0) { // nexus
        std::string tline;
        // are we at the beginning of the file?
        // TODO: add check for interleave and kick out to do a different reader
        // checks for beginning of char by MATRIX
        if (!retstring.empty() && retstring[0] == '#') {
            bool found = false;
            while (getline_safe(stri, tline)) {
                trim_spaces(tline);
                tline = string_to_upper(tline);
                if (tline == "MATRIX") {
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cout << "badly formatted nexus file, missing 'MATRIX' in data/character block" << std::endl;
            }
            retstring = "";
        }
        getline_safe(stri, tline);
        trim_spaces(tline);
        while (tline.empty() || check_comment_line(tline)) {
            if (getline_safe(stri, tline)) {
                trim_spaces(tline);
            } else {
                return false;
            }
        }
        std::vector<std::string> tokens;
        std::string del(" \t");
        tokenize(tline, tokens, del);
        
        // need to check for quoted labels here
        
        if (tokens.size() > 1) {
            for (unsigned int i = 0; i < tokens.size(); i++) {
                trim_spaces(tokens[i]);
            }
            if (tokens[0] == ";") {
                return false;
            } else if (tokens[0][0] == '\'') { // treat ' and " cases separately
                std::string::size_type start = tline.find_first_of('\'');
                std::string::size_type stop  = tline.find_last_of('\'');
                std::string label = tline.substr(start, stop - start + 1);
                std::string seqstr = tline.erase(start, stop - start + 1);
                trim_spaces(seqstr);
                
                seq.set_id(label);
                seq.set_sequence(seqstr);
                return true;
            } else if (tokens[0][0] == '\"') {
                std::string::size_type start = tline.find_first_of('\"');
                std::string::size_type stop  = tline.find_last_of('\"');
                std::string label = tline.substr(start, stop - start + 1);
                std::string seqstr = tline.erase(start, stop - start + 1);
                trim_spaces(seqstr);
                
                seq.set_id(label);
                seq.set_sequence(seqstr);
                return true;
            } else {
                seq.set_id(tokens[0]);
                seq.set_sequence(tokens[1]);
                return true;
            }
        } else {
            return false;
        }
    } else if (ftype == 1) { // phylip
        std::vector<std::string> tokens;
        std::string del(" \t");
        std::string tline;
        // check to see if we are at the beginning of the file
        if (!retstring.empty()) {
            tokenize(retstring, tokens, del);
            if (tokens.size() > 1) {
                trim_spaces(tokens[0]);
                if (is_number(tokens[0])) {
                    getline_safe(stri, tline);
                } else {
                    tline = retstring;
                }
            }
            retstring = "";
        }
        if (tline.empty()) {
            if (!getline_safe(stri, tline)) {
                return false;
            }
            if (tline.empty()) {
                return false;
            }
        }
        tokens.clear();
        tokenize(tline, tokens, del);
        for (unsigned int i = 0; i < tokens.size(); i++) {
            trim_spaces(tokens[i]);
        }
        if (tokens[0].empty()) {
            return false;
        }
        seq.set_id(tokens[0]);
        //TESTING for if there are spaces in the sequences
        // as would be the case for many pxstrec
        if (tokens.size() == 2) {
            seq.set_sequence(tokens[1]);
        } else {
            std::string tse = tokens[1];
            //TODO: look for decimal and add cont char if decimal present
            //seq.add_multistate_char(atoi(tokens[1].c_str()));
            for (unsigned int j=2; j < tokens.size(); j++) {
                tse += " " + tokens[j];
            //seq.add_multistate_char(atoi(tokens[j].c_str()));
            }
            seq.set_sequence(tse);
        }
        return true;
    } else if (ftype == 2) { // fasta
        bool first = true;
        bool going = true;
        std::string curseq;
        while (going) {
            if (first && !retstring.empty()) {
                tline = retstring;
                retstring = "";
            } else {
                if (!getline_safe(stri, tline)) {
                    trim_spaces(curseq);
                    seq.set_sequence(curseq);
                    return false;
                }
            }
            if (tline.substr(0, 1) == ">") {
                if (first) {
                    std::string id_ = tline.substr(1, tline.size()-1);
                    first = false;
                    seq.set_id(id_);
                    curseq = "";
                } else {
                    trim_spaces(curseq);
                    seq.set_sequence(curseq);
                    retstring = tline;
                    return true;
                }
            } else {
                curseq.append(tline);
            }
        }
    } else if (ftype == 3) { // fastq assumes a 33 offset for now
        std::string line1, line2, line3, line4;
        if (!retstring.empty()) {
            line1 = retstring;
            retstring = "";
        } else if (!getline_safe(stri, line1)) {
            return false;
        }
        if (!getline_safe(stri, line2)) {
            return false;
        }
        if (!getline_safe(stri, line3)) {
            return false;
        }
        if (!getline_safe(stri, line4)) {
            return false;
        }
        seq.set_id(line1.substr(1, line1.size()-1));
        seq.set_sequence(line2);
        seq.set_qualstr(line4, 33);
        return true;
    }
    return false;
}


// don't search for MATRIX; if we know it is interleaved, MATRIX has already been read
// need to check for internal comments
// TODO: need to be able to handle quoted strings (specifically, labels that contain spaces)
//       - see function above to see how to do this
std::vector<Sequence> read_interleaved_nexus (std::istream& stri, int num_taxa, int num_char) {
    std::vector<Sequence> seqs;
    std::string tline;
    
    int totcount = 0;
    int loopcount = 0;
    std::string del(" \t");
    while (getline_safe(stri, tline)) {
        trim_spaces(tline);
        if (!tline.empty()) {
            if (check_nexus_comment(tline)) {
                process_nexus_comment(stri, tline);
                continue;
            }
            std::vector<std::string> tokens;
            tokenize(tline, tokens, del);
            if (tokens.size() > 1) {
                Sequence seq;
                for (unsigned int i = 0; i < tokens.size(); i++) {
                    trim_spaces(tokens[i]);
                }
                if (tokens[0] == ";") {
                    std::cout << "Huh?" << std::endl;
                    exit(0);
                } else if (tokens[0][0] == '\'') { // treat ' and " cases separately
                    std::string::size_type start = tline.find_first_of('\'');
                    std::string::size_type stop  = tline.find_last_of('\'');
                    std::string label = tline.substr(start, stop - start + 1);
                    std::string seqstr = tline.erase(start, stop - start + 1);
                    trim_spaces(seqstr);

                    seq.set_id(label);
                    seq.set_sequence(seqstr);
                    
                    if (totcount < num_taxa) {
                        seqs.push_back(seq);
                    } else {
                        // finished with first block. match and concatenate with existing seq
                        // as first pass, assume same ordering (demonic if it isn't)
                        seqs[loopcount].set_sequence(seqs[loopcount].get_sequence() + seq.get_sequence());
                    }
                } else if (tokens[0][0] == '\"') {
                    std::string::size_type start = tline.find_first_of('\"');
                    std::string::size_type stop  = tline.find_last_of('\"');
                    std::string label = tline.substr(start, stop - start + 1);
                    std::string seqstr = tline.erase(start, stop - start + 1);
                    trim_spaces(seqstr);

                    seq.set_id(label);
                    seq.set_sequence(seqstr);
                    
                    if (totcount < num_taxa) {
                        seqs.push_back(seq);
                    } else {
                        // finished with first block. match and concatenate with existing seq
                        // as first pass, assume same ordering (demonic if it isn't)
                        seqs[loopcount].set_sequence(seqs[loopcount].get_sequence() + seq.get_sequence());
                    }
                    
                } else {
                    seq.set_id(tokens[0]);
                    seq.set_sequence(tokens[1]);
                    if (totcount < num_taxa) {
                        seqs.push_back(seq);
                    } else {
                        // finished with first block. match and concatenate with existing seq
                        // as first pass, assume same ordering (demonic if it isn't)
                        seqs[loopcount].set_sequence(seqs[loopcount].get_sequence() + seq.get_sequence());
                    }
                }
            }
            totcount++;
            loopcount++;
            if (loopcount == num_taxa) {
                loopcount = 0; // reset
                // check if we're done
                std::string terp = seqs[num_taxa - 1].get_sequence();
                if ((int)terp.size() == num_char) {
                    break;
                }
            }
        }
    }
    return seqs;
}

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta)
 * TODO: need to add csv
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
*/
// *** this is only used by pxconrates ***
int test_char_filetype_stream(std::istream& stri, std::string& retstring) {
    if (!getline_safe(stri, retstring)) {
        std::cout << "ERROR: end of file too soon" << std::endl;
    }
    int ret = 666; // if you get 666, there is no filetype set
    //NEXUS
    if (retstring[0] == '#') {
        ret = 0;
    } else if (retstring[0] == '>') {
        ret = 2;
    } else {
        std::vector<std::string> tokens;
        std::string del(" \t");
        tokenize(retstring, tokens, del);
        if (tokens.size() > 1) {
            trim_spaces(tokens[0]);
            if (is_number(tokens[0])) {
                ret = 1;
            }
        }
    }
    return ret;
}


/*
 * returns the next string in the getline if there is one
 * TODO: nexus 
 * TODO: decide if this should just be merged with the above reader
*/
// *** this is only used by pxconrates ***
bool read_next_seq_char_from_stream (std::istream& stri, int ftype, std::string& retstring, Sequence& seq) {
    std::string tline;
    if (ftype == 1) { // phylip
        std::vector<std::string> tokens;
        tokens.clear();
        std::string del(" \t");
        std::string tline;
        // check to see if we are at the beginning of the file
        if (!retstring.empty()) {
            tokenize(retstring, tokens, del);
            if (tokens.size() > 1) {
                trim_spaces(tokens[0]);
                if (is_number(tokens[0])) {
                    getline_safe(stri, tline);
                } else {
                    tline = retstring;
                }
            }
            retstring = "";
        }
        if (tline.empty()) {
            if (!getline_safe(stri, tline)) {
                return false;
            }
        }
        tokens.clear();
        tokenize(tline, tokens, del);
        for (unsigned int i = 0; i < tokens.size(); i++) {
            trim_spaces(tokens[i]);
        }
        if (tokens[0].empty()) {
            return false;
        }
        seq.set_id(tokens[0]);
        // split the tokens by spaces
        for (unsigned int i = 1; i < tokens.size(); i++) {
            seq.add_cont_char((double)atof(tokens[i].c_str()));
        }
        return true;
    } else if (ftype == 2) { // fasta
        bool first = true;
        bool going = true;
        std::vector<std::string> tokens;
        std::string del(" \t");
        std::string tline;
        std::string curseq;
        while (going) {
            if (first && !retstring.empty()) {
                tline = retstring;
                retstring = "";
            } else {
                if (!getline_safe(stri, tline)) {
                    tokens.clear();
                    tokenize(curseq, tokens, del);
                    for (unsigned int i = 0; i < tokens.size(); i++) {
                        trim_spaces(tokens[i]);
                    }
                    if (tokens[0].empty()) {
                        return false;
                    }
                    for (unsigned int i = 0; i < tokens.size(); i++) {
                        seq.add_cont_char((double)atof(tokens[i].c_str()));
                    }
                    return false;
                }
            }
            if (tline.substr(0, 1) == ">") {
                if (first) {
                    std::string id_ = tline.substr(1, tline.size()-1);
                    first = false;
                    seq.set_id(id_);
                    curseq = "";
                } else {
                    // split the tokens by spaces
                    tokens.clear();
                    tokenize(curseq, tokens, del);
                    for (unsigned int i = 0; i < tokens.size(); i++) {
                        trim_spaces(tokens[i]);
                    }
                    if (tokens[0].empty()) {
                        return false;
                    }
                    for (unsigned int i = 0; i < tokens.size(); i++) {
                        seq.add_cont_char((double)atof(tokens[i].c_str()));
                    }
                    retstring = tline;
                    return true;
                }
            } else {
                curseq.append(tline);
            }
        }
    }
    return false;
}


// overloaded stream version
void get_nexus_dimensions (std::istream& stri, int& num_taxa, int& num_char, bool& interleave) {
    num_taxa = num_char = 0;
    std::string tline;
    while (getline_safe(stri, tline)) {
        if (!tline.empty()) {
            // check for comments. could be anywhere
            if (check_nexus_comment(tline)) {
                process_nexus_comment(stri, tline);
                continue;
            }
            // convert to uppercase
            tline = string_to_upper(tline);
            std::vector<std::string> searchtokens = tokenize(tline);
            if (searchtokens.empty()) {
                // this will be the case if only whitespace (essentially an empty line)
                continue;
            }
            
            // if we have separate TAXA and CHARACTER blocks, intervening ';' and 'END;' are ignored
            // important thing is to stop when we hit 'MATRIX'
            
            if (searchtokens[0] == "DIMENSIONS") {
            // get rid of '=' and ';'. tokens then easy to deal with.
                std::replace(tline.begin(), tline.end(), '=', ' ');
                std::replace(tline.begin(), tline.end(), ';', ' ');
                searchtokens = tokenize(tline);
                for (unsigned int i = 0; i < searchtokens.size(); i++) {
                    if (searchtokens[i].substr(0, 4) == "NTAX") {
                        i++;
                        num_taxa = stoi(searchtokens[i]);
                    } else if (searchtokens[i].substr(0, 4) == "NCHA") {
                        i++;
                        num_char = stoi(searchtokens[i]);
                    }
                }
            } else if (searchtokens[0] == "FORMAT") {
                std::replace(tline.begin(), tline.end(), '=', ' ');
                std::replace(tline.begin(), tline.end(), ';', ' ');
                searchtokens = tokenize(tline);
                for (unsigned int i = 0; i < searchtokens.size(); i++) {
                    if (searchtokens[i].substr(0, 4) == "INTE") {
                        if (i < (searchtokens.size() - 1)) {
                            i++;
                            if (searchtokens[i] == "YES") {
                                interleave = true;
                            } else if (searchtokens[i] == "NO") {
                                interleave = false;
                            } else {
                                // if yes or no not provided, it is true
                                interleave = true;
                                i--; // backup, since this is a different token
                            }
                        } else {
                            // if `interleave` is the last token, it is true
                            interleave = true;
                        }
                    }
                }
            } else if (searchtokens[0] == "MATRIX") {
                break;
            }
        }
    }
}


// same as above, but grabs datatype and (possibly) 'symbols' (for morphology)
// should remove global to_upper as morphology can be coded arbitrarily
// - this is _low_ priority
// need to deal with files that have a separate taxa block:
// BEGIN TAXA;
//   DIMENSIONS NTAX=1215;
//   TAXLABELS
//     ...
//   ;
// END;
// BEGIN CHARACTERS;
//   DIMENSIONS NCHAR=11395;
//   FORMAT MISSING=? GAP=- DATATYPE=DNA INTERLEAVE=NO;
//   MATRIX
          
void get_nexus_alignment_properties (std::istream& stri, int& num_taxa, int& num_char,
        bool& interleave, std::string& alpha_name, std::string& symbols, char& gap, char& missing) {
    num_taxa = num_char = 0;
    alpha_name = symbols = "";
    // set defaults, in case not explicitly stated
    gap = '-';
    missing = '?';
    
    std::string tline;
    while (getline_safe(stri, tline)) {
        if (!tline.empty()) {
            // check for comments. could be anywhere
            if (check_nexus_comment(tline)) {
                process_nexus_comment(stri, tline);
                continue;
            }
            // convert to uppercase
            tline = string_to_upper(tline);
            std::vector<std::string> searchtokens = tokenize(tline);
            if (searchtokens.empty()) {
                // this will be the case if only whitespace (essentially an empty line)
                continue;
            }
            
            // if we have separate TAXA and CHARACTER blocks, intervening ';' and 'END;' are ignored
            // important thing is to stop when we hit 'MATRIX'
            
            if (searchtokens[0] == "DIMENSIONS") {
            // get rid of '=' and ';'. tokens then easy to deal with.
                std::replace(tline.begin(), tline.end(), '=', ' ');
                std::replace(tline.begin(), tline.end(), ';', ' ');
                searchtokens = tokenize(tline);
                for (unsigned int i = 0; i < searchtokens.size(); i++) {
                    if (searchtokens[i].substr(0, 4) == "NTAX") {
                        i++;
                        num_taxa = stoi(searchtokens[i]);
                    } else if (searchtokens[i].substr(0, 4) == "NCHA") {
                        i++;
                        num_char = stoi(searchtokens[i]);
                    }
                }
            } else if (searchtokens[0] == "FORMAT") {
                std::replace(tline.begin(), tline.end(), '=', ' ');
                std::replace(tline.begin(), tline.end(), ';', ' ');
                searchtokens = tokenize(tline);
                for (unsigned int i = 0; i < searchtokens.size(); i++) {
                    if (searchtokens[i].substr(0, 4) == "INTE") {
                        if (i < (searchtokens.size() - 1)) {
                            i++;
                            if (searchtokens[i] == "YES") {
                                interleave = true;
                            } else if (searchtokens[i] == "NO") {
                                interleave = false;
                            } else {
                                // if yes or no not provided, it is true
                                interleave = true;
                                i--; // backup, since this is a different token
                            }
                        } else {
                            // if `interleave` is the last token, it is true
                            interleave = true;
                        }
                    } else if (searchtokens[i] == "DATATYPE") {
                        i++;
                        // valid Nexus types: STANDARD, DNA, RNA, Nucleotide, Protein, Continuous
                        // phyx types: DNA = 0, AA = 1, BINARY = 2, MULTI = 3, CODON = 4, NA = 5
                        // only considering major ones atm
                        if (searchtokens[i] == "STANDARD") {
                            alpha_name = "MULTI";
                        } else if (searchtokens[i] == "DNA") {
                            alpha_name = "DNA";
                        } else if (searchtokens[i] == "PROTEIN") {
                            alpha_name = "AA";
                        } else {
                            std::cout << "Datatype '" << searchtokens[i] << "' not supported" << std::endl;
                        }
                    } else if (searchtokens[i] == "GAP") {
                        i++;
                        gap = searchtokens[i][0];
                    } else if (searchtokens[i] == "MISSING") {
                        i++;
                        missing = searchtokens[i][0];
                    } else if (searchtokens[i] == "SYMBOLS") {
                        // morphology data
                        // can take form SYMBOLS="012" or (annoyingly) SYMBOLS="0 1 2"
                        // just need to make sure both " are captured, and we should be cool
                        bool done = false;
                        bool first = true;
                        std::string terp;
                        while (!done) {
                            i++;
                            if (first) {
                                // check that we are in the right spot
                                terp = searchtokens[i];
                                if (terp.front() == '"') {
                                    // check if symbols contain no gaps
                                    if (terp.back() == '"') {
                                        //std::cout << "symbols are (rationally) contiguous!" << std::endl;
                                        terp.erase (std::remove (terp.begin(), terp.end(), '"'), terp.end());
                                        symbols = terp;
                                        done = true;
                                    }
                                }
                                first = false;
                            } else {
                                terp += searchtokens[i];
                                if (searchtokens[i].back() == '"') {
                                    //std::cout << "found the end of symbols!" << std::endl;
                                    terp.erase (std::remove (terp.begin(), terp.end(), '"'), terp.end());
                                    symbols = terp;
                                    done = true;
                                }
                            }
                        }
                        //std::cout << "Captured symbols: " << symbols << std::endl;
                    }
                }
            } else if (searchtokens[0] == "MATRIX") {
                if (alpha_name == "MULTI" && symbols.empty()) {
                    symbols = "01"; // this is default Nexus alphabet for 'standard' data
                    alpha_name = "BINARY";
                }
                break;
            }
        }
    }
}


void get_phylip_dimensions (const std::string& head, int& num_taxa, int& num_char) {
    std::vector<std::string> header = tokenize(head);
    num_taxa = stoi(header[0]);
    num_char = stoi(header[1]);
}


bool is_complicated_phylip (std::istream& pios, const int& num_char) {
    bool complicated = false;
    std::string peek = peek_line(pios);
    std::vector<std::string> tokens = tokenize(peek);
    
    if (tokens.size() != 2) {
        // if only 1, no space between label and sequence
        // if > 2, internal spaces in the sequence itself
        complicated = true;
    } else {
        // if it seems good, check that length of seq equals stated value
        // if not, we are dealing with 1) interleaved or 2) multi-line
        if ((int)tokens[1].size() != num_char) {
            complicated = true;
        }
    }
    return complicated;
}


/* figure out non-simple phylip alignments formats (i.e. not contiguous):
1) multiline data (like exported from paup). e.g.:
    5 20
    TaxonA  AAATTTCCCTGTCCC
            TTTAA
    TaxonB  GCTCGAGGGGCCCCA
            AGACC
    TaxonC  ACGCTCCCCCTTAAA
            AATGA
    TaxonD  TCCTTGTTCAACTCC
            GGTGG
    TaxonE  TTACTATTCCCCCCC
            GCCGG
2) interleaved data. e.g.:
    5 20
    TaxonA    AAATTTCCCTGTCCC
    TaxonB    GCTCGAGGGGCCCCA
    TaxonC    ACGCTCCCCCTTAAA
    TaxonD    TCCTTGTTCAACTCC
    TaxonE    TTACTATTCCCCCCC
                    <- blank line may or may not be present
    TTTAA
    AGACC
    AATGA
    GGTGG
    GCCGG
3) sequence includes spaces. e.g.:
    5 20
    TaxonA  AAATT TCCCT GTCCC TTTAA
    TaxonB  GCTCG AGGGG CCCCA AGACC
    TaxonC  ACGCT CCCCC TTAAA AATGA
    TaxonD  TCCTT GTTCA ACTCC GGTGG
    TaxonE  TTACT ATTCC CCCCC GCCGG
- header has already been processed (num_taxa, num_char)
- all bools have been initialized to false above
- for spaces-only, we will know the format after reading the first line.
- atm we are _not_ considering 10 character limit labels
  - i.e., such tht there may not be a space between the label and sequence
  - this includes not allowing internal spaces in taxon labels
- also assuming that there is not a combination of multiline _and_ interleaved
*/
void get_phylip_format (std::istream& pios, const unsigned int& num_taxa, const unsigned int& num_char,
        bool& interleaved, bool& spaces, bool& multiline) {
    // store current position of the stream so we can rewind after determining format
    std::streampos spt = pios.tellg();
    
    std::string line;
    std::string name;
    std::string seq;
    int num_elem = 0; // count chunks when spaces present
    std::vector<std::string> tokens;
    
    // debug stuff
    //bool done = false;
    //int lcnt = 0; // line counter for multiline seqs
    
    // first line
    getline_safe(pios, line);
    tokens = tokenize(line);
    name = tokens[0];
    num_elem = (int)tokens.size();
    if (num_elem == 2) {
        // could be: 1) simple, 2) multiline, or interleaved
        seq = tokens[1];
        if (seq.length() == num_char) {
            // return to the original position in the stream
            pios.seekg(spt, std::ios_base::beg);
            return; // simple single line //
        } else {
            // read in next line
            getline_safe(pios, line);
            tokens = tokenize(line);
            if (tokens.size() == 1) {
                // simple multiline format
                multiline = true;
                // keep reading to check
                /*
                done = false;
                while (!done) {
                    seq += tokens[0];
                    if (seq.length() == num_char) {
                        done = true;
                    } else {
                        getline_safe(pios, line);
                        tokens = tokenize(line);
                    }
                }
                std::cout << "Successfully read in multiline sequence ;)" << std::endl;
                std::cout << seq << std::endl;
                */
                // return to the original position in the stream
                pios.seekg(spt, std::ios_base::beg);
                return; // multiline //
            } else if (tokens.size() == 2) {
                // likely interleaved format
                interleaved = true;
                // keep reading to check. seq continues every num_taxa lines
                /*
                done = false;
                int tcnt = 2; // we've read in two already
                lcnt = 1;
                while (!done) {
                    getline_safe(pios, line);
                    tcnt++;
                    //std::cout << "tcnt = " << tcnt << std::endl;
                    if (tcnt == num_taxa) {
                        //std::cout << "skipping line: " << line << std::endl;
                        tcnt = 0;
                    } else if (tcnt == 1) {
                        lcnt++;
                        if (line.empty()) {
                            // there may be an empty line between blocks
                            //std::cout << "skipping an empty line." << std::endl;
                            getline_safe(pios, line);
                        }
                        tokens = tokenize(line);
                        if (tokens.size() > 1) {
                            std::cout << "Um, I'm not sure what is happening..." << std::endl;
                        } else {
                            seq += tokens[0];
                            if (seq.length() == num_char) {
                                done = true;
                            }
                        }
                    } else {
                        //std::cout << "skipping line: " << line << std::endl;
                    }
                }
                std::cout << "Successfully read in interleaved sequence ("
                        << lcnt << " lines) ;)" << std::endl;
                std::cout << seq << std::endl;
                */
                // return to the original position in the stream
                pios.seekg(spt, std::ios_base::beg);
                return; // interleaved //
            } else {
                std::cerr << "Error: unrecognized phylip format. Exiting." << std::endl; 
                exit(1);
            }
        }
    } else if (num_elem > 2) {
        // dealing with spaces
        spaces = true;
        for (unsigned int i = 1; i < tokens.size(); i++) {
            seq += tokens[i];
        }
        // check if that is all that is going on
        if (seq.length() == num_char) {
            //std::cout << "Successfully read in spaces-only sequence ;)" << std::endl;
            //std::cout << seq << std::endl;
            // return to the original position in the stream
            pios.seekg(spt, std::ios_base::beg);
            return; // just spaces //
        } else {
            // now, determine if multiline or interleaved
            // if multiline, 2nd line will have (at maximum) 1 fewer element
            // if interleaved, 2nd line will have the same number of elements
            getline_safe(pios, line);
            tokens = tokenize(line);
            if ((int)tokens.size() == num_elem) {
                interleaved = true;
                // keep reading to check. seq continues every num_taxa lines
                /*
                done = false;
                int tcnt = 2; // we've read in two already
                lcnt = 1;
                while (!done) {
                    getline_safe(pios, line);
                    tcnt++;
                    //std::cout << "tcnt = " << tcnt << std::endl;
                    if (tcnt == num_taxa) {
                        //std::cout << "skipping line: " << line << std::endl;
                        tcnt = 0;
                    } else if (tcnt == 1) {
                        lcnt++;
                        if (line.empty()) {
                            // there may be an empty line between blocks
                            //std::cout << "skipping an empty line." << std::endl;
                            getline_safe(pios, line);
                        }
                        tokens = tokenize(line);
                        for (unsigned int i = 0; i < tokens.size(); i++) {
                            seq += tokens[i];
                        }
                        if (seq.length() == num_char) {
                            done = true;
                        }
                    } else {
                        //std::cout << "skipping line: " << line << std::endl;
                    }
                }
                std::cout << "Successfully read in spaces+interleaved sequence ("
                        << lcnt << " lines) ;)" << std::endl;
                std::cout << seq << std::endl;
                */
                // return to the original position in the stream
                pios.seekg(spt, std::ios_base::beg);
                return; // spaces & interleaved //
            } else {
                multiline = true;
                // check to make sure
                /*
                lcnt = 1;
                done = false;
                while (!done) {
                    for (unsigned int i = 0; i < tokens.size(); i++) {
                        seq += tokens[i];
                        lcnt++;
                    }
                    if (seq.length() == num_char) {
                        done = true;
                    }
                }
                std::cout << "Successfully read in spaces+multiline sequence ("
                        << lcnt << " lines) ;)" << std::endl;
                std::cout << seq << std::endl;
                */
                // return to the original position in the stream
                pios.seekg(spt, std::ios_base::beg);
                return; // spaces & multiline //
            }
        }
    } else {
        std::cerr << "Error: unrecognized phylip format. Exiting." << std::endl; 
        exit(1);
    }
}


// read in various phylip formats
// do you _see_ how needlessly complicated this is JF?!?
std::vector<Sequence> read_phylip (std::istream& pios, const int& num_taxa, const int& num_char) {
    
    std::vector<Sequence> seqs;
    Sequence seq;
    std::string line;
    bool spaces = false;
    bool multiline = false;
    bool interleaved = false;
    bool first = true;
    int lineCount = 0;
    //int charCount = 0;
    std::vector<std::string> tokens;
    std::vector<std::string> names;
    std::string residues;
    std::string name;
    
    get_phylip_format(pios, num_taxa, num_char, interleaved, spaces, multiline);
    
    //std::cout << "interleaved = " << interleaved << "; spaces = "
    //        << spaces << "; multiline = " << multiline<< std::endl;
    if ((spaces + multiline + interleaved) == 0) {
        // simple format. standard read
        //std::cout << std::endl << "processing simple contiguous sequencess" << std::endl;
        while (read_next_seq_from_stream(pios, 1, line, seq)) {
            seqs.push_back(seq);
        }
    } else if (multiline) {
        if (spaces) {
            //std::cout << std::endl << "processing multiline sequences with spaces" << std::endl;
            while (getline_safe(pios, line)) {
                if (line.empty()) {
                    continue;
                }
                tokens = tokenize(line);
                if (first) {
                    lineCount++;
                    name = tokens[0];
                    for (unsigned int i = 1; i < tokens.size(); i++) {
                        residues += tokens[i];
                    }
                    first = false;
                } else {
                    for (unsigned int i = 0; i < tokens.size(); i++) {
                        residues += tokens[i];
                    }
                }
                if ((int)residues.size() == num_char) {
                    seq.set_id(name);
                    seq.set_sequence(residues);
                    seqs.push_back(seq);
                    residues = "";
                    name = "";
                    first = true;
                }
            }
        } else {
            //std::cout << std::endl << "processing multiline sequences" << std::endl;
            while (getline_safe(pios, line)) {
                if (line.empty()) {
                    continue;
                }
                tokens = tokenize(line);
                if (first) {
                    lineCount++;
                    name = tokens[0];
                    residues = tokens[1];
                    first = false;
                } else {
                    residues += tokens[0];
                }
                if ((int)residues.size() == num_char) {
                    seq.set_id(name);
                    seq.set_sequence(residues);
                    seqs.push_back(seq);
                    residues = "";
                    name = "";
                    first = true;
                }
            }
        }
    } else if (interleaved) {
        if (spaces) {
            //std::cout << std::endl << "processing interleaved sequences with spaces" << std::endl;
            int taxcnt = 0;
            lineCount = 0;
            while (getline_safe(pios, line)) {
                if (line.empty()) {
                    continue;
                }
                tokens = tokenize(line);
                if (lineCount < num_taxa) {
                    name = tokens[0];
                    for (unsigned int i = 1; i < tokens.size(); i++) {
                        residues += tokens[i];
                    }
                    seq.set_id(name);
                    seq.set_sequence(residues);
                    seqs.push_back(seq);
                    residues = "";
                    name = "";
                } else {
                    for (unsigned int i = 0; i < tokens.size(); i++) {
                        residues += tokens[i];
                    }
                    seqs[taxcnt].set_sequence(seqs[taxcnt].get_sequence() + residues);
                    taxcnt++;
                    if (taxcnt == num_taxa) {
                        taxcnt = 0;
                    }
                    residues = "";
                }
                lineCount++;
            }
        } else {
            //std::cout << std::endl << "processing interleaved sequences" << std::endl;
            int taxcnt = 0;
            lineCount = 0;
            while (getline_safe(pios, line)) {
                if (line.empty()) {
                    continue;
                }
                tokens = tokenize(line);
                if (lineCount < num_taxa) {
                    name = tokens[0];
                    residues = tokens[1];
                    seq.set_id(name);
                    seq.set_sequence(residues);
                    seqs.push_back(seq);
                } else {
                    seqs[taxcnt].set_sequence(seqs[taxcnt].get_sequence() + tokens[0]);
                    taxcnt++;
                    if (taxcnt == num_taxa) {
                        taxcnt = 0;
                    }
                }
                lineCount++;
            }
        }
        // check seq lengths are all good
        bool good = true;
        int bad_idx = 0;
        int bad_len = 0;
        for (unsigned int i = 0; i < seqs.size(); i++) {
            residues = seqs[i].get_sequence();
            if ((int)residues.size() != num_char) {
                bad_idx = i;
                bad_len = (int)residues.size();
                good = false;
                break; // break after first failure
            }
        }
        if (!good) {
            std::cerr << "Error: bad phylip file. Taxon '" << seqs[bad_idx].get_id()
                << "' had " << bad_len << " characters, but " << num_char
                << " were expected. Exiting." << std::endl;
            exit(1);
        }
    } else if (spaces) {
        //std::cout << std::endl << "processing sequences with spaces" << std::endl;
        while (getline_safe(pios, line)) {
            if (line.empty()) {
                continue;
            }
            tokens = tokenize(line);
            name = tokens[0];
            for (unsigned int i = 1; i < tokens.size(); i++) {
                residues += tokens[i];
            }
            if ((int)residues.size() != num_char) {
                std::cerr << "Error: bad phylip file. Taxon '" << name
                    << "' had " << residues.size() << " characters, but " << num_char
                    << " were expected. Exiting." << std::endl;
                    exit(1);
            }
            seq.set_id(name);
            seq.set_sequence(residues);
            seqs.push_back(seq);
            residues = "";
        }
    } else {
        std::cerr << "Error: unrecognized phylip format. Exiting." << std::endl; 
        exit(1);
    }
    
    //std::cout << "Read in " << seqs.size() << " sequences." << std::endl;
    if ((int)seqs.size() != num_taxa) {
        std::cerr << "Error: bad phylip file. Read in " << seqs.size()
                << " sequences, but expected " << num_taxa << ". Exiting." << std::endl;
        exit(1);
    }
    return seqs;
}


// general reader that is used by a bunch of stuff (i.e., avoid code duplication)
// returns vector of sequences and (by reference) alphabet name
// used by programs that need the entire alignment in memory at once
// all other properties are determined elsewhere
std::vector<Sequence> ingest_alignment (std::istream* pios, std::string& alphaName) {
    std::vector<Sequence> seqs;
    Sequence seq;
    std::string retstring;
    alphaName = "";
    int ft = test_seq_filetype_stream(*pios, retstring);
    int file_num_taxa = 0; // num_taxa declared in the file itself
    int file_num_char = 0; // likewise for num_char
    std::string file_type = get_filetype_string(ft);
    
    if (file_type == "nexus") {
        bool is_interleaved = false;
        // a bunch of required variables, not used here:
        char gap, missing;
        std::string symbols;
        get_nexus_alignment_properties(*pios, file_num_taxa, file_num_char,
                is_interleaved, alphaName, symbols, gap, missing);
        //std::cout << "alphaName = " << alphaName << std::endl;
        retstring = ""; // have to set so seq_reader knows we are mid-file
        if (!is_interleaved) {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                seqs.push_back(seq);
            }
        } else {
            seqs = read_interleaved_nexus(*pios, file_num_taxa, file_num_char);
        }
    } else {
        bool complicated_phylip = false;
        // check if we are dealing with a complicated phylip format
        if (file_type == "phylip") {
            get_phylip_dimensions(retstring, file_num_taxa, file_num_char);
            complicated_phylip = is_complicated_phylip(*pios, file_num_char);
        }
        if (complicated_phylip) {
            seqs = read_phylip(*pios, file_num_taxa, file_num_char);
            if (alphaName.empty()) {
                alphaName = seqs[0].get_alpha_name();
            }
        } else {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                seqs.push_back(seq);
                if (alphaName.empty()) {
                    alphaName = seq.get_alpha_name();
                }
            }
            if (ft == 2) { // fasta has an trailing one
                seqs.push_back(seq);
            }
        }
    }
    
    // some simple error-checking
    if (file_num_taxa != 0) {
        if (file_num_taxa != (int)seqs.size()) {
            std::cerr << "Error: number of taxa declared in the file ("
                << ") does not match the number read (" << seqs.size()
                << "). Exiting." << std::endl;
            exit(1);
        }
    }
    if (file_num_char != 0) {
        // if num_char comes from a file, _must_ be aligned
        bool aligned = is_aligned(seqs);
        if (!aligned) {
            std::cerr << "Error: sequences are not aligned. Exiting." << std::endl;
            exit(0);
        }
    }
    
    return seqs;
}
