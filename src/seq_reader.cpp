#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <map>

#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"


// for printing purposes
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
    if (!getline(stri, retstring)) {
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
 * TODO: nexus interleaved is not going to work here
 * TODO: skip Nexus comment lines
 */
bool read_next_seq_from_stream (std::istream & stri, int ftype, std::string& retstring, Sequence& seq) {
    std::string tline;
    if (ftype == 0) { // nexus
        std::string tline;
        //are we at the beginning of the file?
        //TODO: add check for interleave and kick out to do a different reader
        //checks for beginning of char by MATRIX
        if (retstring.size() > 0 && retstring[0] == '#') {
            bool found = false;
            while (getline(stri, tline)) {
                trim_spaces(tline);
                tline = string_to_upper(tline);
                if (tline.compare("MATRIX") == 0) {
                    found = true;
                    break;
                }
            }
            if (found == false) {
                std::cout << "badly formatted nexus file, missing 'MATRIX' in data/character block" << std::endl;
            }
            retstring = "";
        }
        getline(stri, tline);
        trim_spaces(tline);
        while (tline.size() == 0  || check_comment_line(tline)) {
            if (getline(stri,tline)) {
                trim_spaces(tline);
            } else {
                return false;
            }
        }
        std::vector<std::string> tokens;
        std::string del(" \t");
        tokenize(tline, tokens, del);
        if (tokens.size() > 1) {
            for (unsigned int i=0; i < tokens.size(); i++) {
                trim_spaces(tokens[i]);
            }
            if (tokens[0].compare(";") == 0) {
                return false;
            } else {
                seq.set_id(tokens[0]);
                seq.set_sequence(tokens[1]);
                return true;
            }
        } else {
            return false;
        }
    } else if (ftype == 1) { // phylip. TODO: this crashes if has a trailing empty line (JWB)
        std::vector<std::string> tokens;
        std::string del(" \t");
        std::string tline;
        // check to see if we are at the beginning of the file
        if (retstring.size() > 0) {
            tokenize(retstring, tokens, del);
            if (tokens.size() > 1) {
                trim_spaces(tokens[0]);
                if (is_number(tokens[0])) {
                    getline(stri, tline);
                } else {
                    tline = retstring;
                }
            }
            retstring = "";
        }
        if (tline.size() == 0) {
            if (!getline(stri, tline)) {
                return false;
            }
            if (tline.size() == 0) {
                return false;
            }
        }
        tokens.clear();
        tokenize(tline, tokens, del);
        for (unsigned int i=0; i < tokens.size(); i++) {
            trim_spaces(tokens[i]);
        }
        if (tokens[0].size() == 0) {
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
        std::string curseq = "";
        while (going) {
            if (first == true && retstring.size() > 0) {
                tline = retstring;
                retstring = "";
            } else {
                if (!getline(stri, tline)) {
                    trim_spaces(curseq);
                    seq.set_sequence(curseq);
                    return false;
                }
            }
            if (tline.substr(0,1) == ">") {
                if (first == true) {
                    std::string id_ = tline.substr(1,tline.size()-1);
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
    } else if (ftype == 3) {//fastq assumes a 33 offset for now
        std::string line1,line2,line3,line4;
        if (retstring.size() > 0) {
            line1 = retstring;
            retstring = "";
        } else if (!getline(stri,line1)) {
            return false;
        }
        if (!getline(stri,line2)) {
            return false;
        }
        if (!getline(stri,line3)) {
            return false;
        }
        if (!getline(stri,line4)) {
            return false;
        }
        seq.set_id(line1.substr(1,line1.size()-1));
        seq.set_sequence(line2);
        seq.set_qualstr(line4,33);
        return true;
    }
    return false;
}


// interleaved data do not work with the stream philosophy
// by using this function, the file has already been checked, so we know ntax and nchar
// prolly get rid of this in favour of the stream-based one
std::vector<Sequence> read_interleaved_nexus_file (std::string filen, int ntax, int nchar) {
    std::vector <Sequence> seqs;
    std::string tline;
    std::ifstream infile(filen.c_str());
    //bool done = false; // not used
    
    // first, get us to the MATRIX line i.e., right before the sequences start
    bool found = false;
    while (getline(infile, tline)) {
        trim_spaces(tline);
        tline = string_to_upper(tline);
        if (tline.compare("MATRIX") == 0) {
            found = true;
            break;
        }
    }
    if (found == false) {
        std::cout << "badly formatted nexus file: missing 'MATRIX' in data/character block. Exiting." << std::endl;
        exit(1);
    }
    
    int totcount = 0;
    int loopcount = 0;
    std::string del(" \t");
    while (getline(infile, tline)) {
        trim_spaces(tline);
        if (tline.size() != 0) {
            std::vector<std::string> tokens;
            tokenize(tline, tokens, del);
            if (tokens.size() > 1) {
                Sequence seq;
                for (unsigned int i=0; i < tokens.size(); i++) {
                    trim_spaces(tokens[i]);
                }
                if (tokens[0].compare(";") == 0) {
                    std::cout << "Huh?" << std::endl;
                    exit(0);
                } else {
                    seq.set_id(tokens[0]);
                    seq.set_sequence(tokens[1]);
                    if (totcount < ntax) {
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
            if (loopcount == ntax) {
                loopcount = 0; // reset
                // check if we're done
                std::string terp = seqs[ntax - 1].get_sequence();
                if ((int)terp.size() == nchar) {
                    break;
                }
            }
        }
    }
    infile.close();
    //std::cout << "Seqs has " << seqs.size() << " taxa and "
    //        << seqs[0].get_sequence().size() << " sites." << std::endl;
    return seqs;
}


// don't search for MATRIX; if we know it is interleaved, MATRIX has already been read
std::vector<Sequence> read_interleaved_nexus (std::istream& stri, int ntax, int nchar) {
    std::vector<Sequence> seqs;
    std::string tline;
    //bool done = false; // not used
    
    int totcount = 0;
    int loopcount = 0;
    std::string del(" \t");
    while (getline(stri, tline)) {
        trim_spaces(tline);
        if (tline.size() != 0) {
            std::vector<std::string> tokens;
            tokenize(tline, tokens, del);
            if (tokens.size() > 1) {
                Sequence seq;
                for (unsigned int i=0; i < tokens.size(); i++) {
                    trim_spaces(tokens[i]);
                }
                if (tokens[0].compare(";") == 0) {
                    std::cout << "Huh?" << std::endl;
                    exit(0);
                } else {
                    seq.set_id(tokens[0]);
                    seq.set_sequence(tokens[1]);
                    if (totcount < ntax) {
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
            if (loopcount == ntax) {
                loopcount = 0; // reset
                // check if we're done
                std::string terp = seqs[ntax - 1].get_sequence();
                if ((int)terp.size() == nchar) {
                    break;
                }
            }
        }
    }
    //std::cout << "Seqs has " << seqs.size() << " taxa and "
    //        << seqs[0].get_sequence().size() << " sites." << std::endl;
    return seqs;
}

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta)
 * TODO: need to add csv
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
 */
int test_char_filetype_stream(std::istream& stri, std::string& retstring) {
    if (!getline(stri, retstring)) {
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
        tokenize(retstring,tokens,del);
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
bool read_next_seq_char_from_stream (std::istream& stri, int ftype, std::string& retstring, Sequence& seq) {
    std::string tline;
    if (ftype == 1) { // phylip
        std::vector<std::string> tokens;
        tokens.clear();
        std::string del(" \t");
        std::string tline;
        //check to see if we are at the beginning of the file
        if (retstring.size() > 0) {
            tokenize(retstring,tokens,del);
            if (tokens.size() > 1) {
                trim_spaces(tokens[0]);
                if (is_number(tokens[0])) {
                    getline(stri,tline);
                } else {
                    tline = retstring;
                }
            }
            retstring = "";
        }
        if (tline.size() == 0) {
            if (!getline(stri,tline)) {
                return false;
            }
        }
        tokens.clear();
        tokenize(tline,tokens,del);
        for (unsigned int i=0; i < tokens.size(); i++) {
            trim_spaces(tokens[i]);
        }
        if (tokens[0].size() == 0) {
            return false;
        }
        seq.set_id(tokens[0]);
        //split the tokens by spaces
        for (unsigned int i=1 ; i < tokens.size(); i++) {
            seq.add_cont_char((double)atof(tokens[i].c_str()));
        }
        return true;
    } else if (ftype == 2) { // fasta
        bool first = true;
        bool going = true;
        std::vector<std::string> tokens;
        std::string del(" \t");
        std::string tline;
        std::string curseq = "";
        while (going) {
            if (first == true && retstring.size() > 0) {
                tline = retstring;
                retstring = "";
            } else {
                if (!getline(stri, tline)) {
                    tokens.clear();
                    tokenize(curseq,tokens,del);
                    for (unsigned int i=0; i < tokens.size(); i++) {
                        trim_spaces(tokens[i]);
                    }
                    if (tokens[0].size() == 0) {
                        return false;
                    }
                    for (unsigned int i=0 ; i < tokens.size(); i++) {
                        seq.add_cont_char((double)atof(tokens[i].c_str()));
                    }
                    return false;
                }
            }
            if (tline.substr(0,1) == ">") {
                if (first == true) {
                    std::string id_ = tline.substr(1,tline.size()-1);
                    first = false;
                    seq.set_id(id_);
                    curseq = "";
                } else {
                    //split the tokens by spaces
                    tokens.clear();
                    tokenize(curseq,tokens,del);
                    for (unsigned int i=0; i < tokens.size(); i++) {
                        trim_spaces(tokens[i]);
                    }
                    if (tokens[0].size() == 0) {
                        return false;
                    }
                    for (unsigned int i =0 ; i < tokens.size(); i++) {
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

// file-version of above. used by concatenator
void get_nexus_dimensions_file (std::string& filen, int& numTaxa, int& numChar, bool& interleave) {
    numTaxa = numChar = 0;
    std::string tline;
    std::string temp;
    std::ifstream infile(filen.c_str());
    while (getline(infile, tline)) {
        if (!tline.empty()) {
            // convert to uppercase
            tline = string_to_upper(tline);
            std::vector<std::string> searchtokens = tokenize(tline);
            if (searchtokens[0] == "DIMENSIONS") {
            // get rid of '=' and ';'. tokens then easy to deal with.
                std::replace(tline.begin(), tline.end(), '=', ' ');
                std::replace(tline.begin(), tline.end(), ';', ' ');
                searchtokens = tokenize(tline);
                for (unsigned int i = 0; i < searchtokens.size(); i++) {
                    if (searchtokens[i].substr(0, 4) == "NTAX") {
                        i++;
                        numTaxa = stoi(searchtokens[i]);
                    } else if (searchtokens[i].substr(0, 4) == "NCHA") {
                        i++;
                        numChar = stoi(searchtokens[i]);
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
    infile.close();
}


// overloaded stream version
void get_nexus_dimensions (std::istream& stri, int& numTaxa, int& numChar, bool& interleave) {
    numTaxa = numChar = 0;
    std::string tline;
    //string temp;
    while (getline(stri, tline)) {
        if (!tline.empty()) {
            // convert to uppercase
            tline = string_to_upper(tline);
            std::vector<std::string> searchtokens = tokenize(tline);
            if (searchtokens[0] == "DIMENSIONS") {
            // get rid of '=' and ';'. tokens then easy to deal with.
                std::replace(tline.begin(), tline.end(), '=', ' ');
                std::replace(tline.begin(), tline.end(), ';', ' ');
                searchtokens = tokenize(tline);
                for (unsigned int i = 0; i < searchtokens.size(); i++) {
                    if (searchtokens[i].substr(0, 4) == "NTAX") {
                        i++;
                        numTaxa = stoi(searchtokens[i]);
                    } else if (searchtokens[i].substr(0, 4) == "NCHA") {
                        i++;
                        numChar = stoi(searchtokens[i]);
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
void get_nexus_alignment_properties (std::istream& stri, int& numTaxa, int& numChar,
        bool& interleave, std::string& alpha_name, std::string& symbols, char& gap, char& missing) {
    numTaxa = numChar = 0;
    alpha_name = symbols = "";
    // set defaults, in case not explicitly stated
    gap = '-';
    missing = '?';
    
    std::string tline;
    //string temp;
    while (getline(stri, tline)) {
        if (!tline.empty()) {
            // convert to uppercase
            tline = string_to_upper(tline);
            std::vector<std::string> searchtokens = tokenize(tline);
            if (searchtokens[0] == "DIMENSIONS") {
            // get rid of '=' and ';'. tokens then easy to deal with.
                std::replace(tline.begin(), tline.end(), '=', ' ');
                std::replace(tline.begin(), tline.end(), ';', ' ');
                searchtokens = tokenize(tline);
                for (unsigned int i = 0; i < searchtokens.size(); i++) {
                    if (searchtokens[i].substr(0, 4) == "NTAX") {
                        i++;
                        numTaxa = stoi(searchtokens[i]);
                    } else if (searchtokens[i].substr(0, 4) == "NCHA") {
                        i++;
                        numChar = stoi(searchtokens[i]);
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
                        std::string terp = "";
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


// **************************** //
// *** Deprecated functions *** //
// **************************** //

// Leftover file-centric functions from before working with streams

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta), @ (fastq)
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
 */

int test_seq_filetype (std::string filen) {
    std::string tline;
    std::ifstream infile(filen.c_str());
    int ret = 666; // if you get 666, there is no filetype set
    while (getline(infile,tline)) {
        if (tline.size() < 1) {
            continue;
        }
        if (tline[0] == '#') {
            ret = 0;
            break;
        }
        std::vector<std::string> tokens;
        std::string del(" \t");
        tokenize(tline,tokens,del);
        if (tokens.size() > 1) {
            trim_spaces(tokens[0]);
            if (is_number(tokens[0])) {
                ret = 1;
                break;
            }
        }
        if (tline[0] == '>') {
            ret = 2;
            break;
        }
        if (tline[0] == '@') {
            ret = 3;
            break;
        }
        break;
    }
    infile.close();
    return ret;
}


/*
bool read_phylip_file(string filen, vector <Sequence>& seqs) {
    string tline;
    ifstream infile(filen.c_str());
    bool first = true;
    while (getline(infile, tline)) {
        vector<string> searchtokens;
        tokenize(tline, searchtokens, "\t    ");
        for (unsigned int j=0; j < searchtokens.size(); j++) {
            trim_spaces(searchtokens[j]);
        }
        if (first == true) { //reading past the first line
            first = false;
            try{
                //int nseqs = atoi(searchtokens[0].c_str());
                //int nsites = atoi(searchtokens[1].c_str());
            }catch( char * str ) {
                return false; //not a phylip
            }
            continue;
        }
        if (searchtokens.size() < 2) {
            continue;
        }
        Sequence a = Sequence(searchtokens[0],searchtokens[1],true);
        seqs.push_back(a);
    }
    infile.close();
    return true;
}
*/

// return false if not a fasta
/*
bool read_fasta_file(string filen, vector <Sequence>& seqs) {
    string tline;
    ifstream infile(filen.c_str());
    bool first = true;
    Sequence cur;
    string curseq;
    while (getline(infile, tline)) {
        trim_spaces(tline);
        if (tline.size() < 1) {
            continue;
        }
        if (tline.substr(0,1) == ">") {
            string id_ = tline.substr(1,tline.size()-1);
            if (first == true) {
                first = false;
                cur = Sequence();
                cur.set_id(id_);
                curseq = "";
            } else {
                cur.set_sequence(curseq);
                seqs.push_back(cur);
                cur = Sequence();
                cur.set_id(id_);
                curseq = "";
            }
        } else {
            curseq += tline;
        }
    }
    cur.set_sequence(curseq);
    seqs.push_back(cur);
    infile.close();
    return true;
}
*/

/*
bool read_phylip_file_strec(string filen, vector <Sequence>& seqs) {
    string tline;
    ifstream infile(filen.c_str());
    bool first = true;
    while (getline(infile, tline)) {
        vector<string> searchtokens;
        tokenize(tline, searchtokens, "\t");
        for (unsigned int j=0; j < searchtokens.size(); j++) {
            trim_spaces(searchtokens[j]);
        }
        if (first == true) { //reading past the first line
            first = false;
            try{
                //int nseqs = atoi(searchtokens[0].c_str());
                //int nsites = atoi(searchtokens[1].c_str());
            }catch( char * str ) {
                return false; //not a phylip
            }
            continue;
        }
        if (searchtokens.size() < 2) {
            continue;
        }
        std::cout << searchtokens[0] << " " << searchtokens[1] << std::endl;
        Sequence a = Sequence(searchtokens[0],searchtokens[1],true);
        seqs.push_back(a);
    }
    infile.close();
    return true;
}
*/

//TODO: INCOMPLETE
/*
bool read_nexus_seqs_file(string filen, vector <Sequence>& seqs) {
    string tline;
    ifstream infile(filen.c_str());
    bool first = true;
    while (getline(infile, tline)) {
        vector <string> searchtokens;
        tokenize(tline, searchtokens, "    ");
        for (unsigned int j=0; j < searchtokens.size(); j++) {
            trim_spaces(searchtokens[j]);
        }
        if (first == true) { //reading past the first line
            first = false;
            try {
                //int nseqs = atoi(searchtokens[0].c_str());
                //int nsites = atoi(searchtokens[1].c_str());
            } catch( char * str ) {
                return false; //not a phylip
            }
            continue;
        }
        if (searchtokens.size() < 2) {
            continue;
        }
        Sequence a = Sequence(searchtokens[0], searchtokens[1], true);
        seqs.push_back(a);
    }
    infile.close();
    return true;
}
*/
