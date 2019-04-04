
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <map>

using namespace std;

#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"

// for printing purposes
map <int, string> filetype_map = {
   {0, "nexus"},
   {1, "phylip"},
   {2, "fasta"},
   {3, "fastq"}
};

string get_filetype_string (int const& ft) {
    string ftype = filetype_map[ft];
    return ftype;
}

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta), @ (fastq)
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
 */
int test_seq_filetype_stream(istream & stri, string & retstring) {
    if (!getline(stri, retstring)) {
        cout << "ERROR: end of file too soon" << endl;
    }
    int ret = 666; // if you get 666, there is no filetype set
    if (retstring[0] == '#') {
        ret = 0;
    } else if (retstring[0] == '>') {
        ret = 2;
    } else if (retstring[0] == '@') {
        ret = 3;
    } else {
        vector <string> tokens;
        string del(" \t");
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
bool read_next_seq_from_stream(istream & stri, int ftype, string & retstring, Sequence & seq) {
    string tline;
    if (ftype == 0) { // nexus
        string tline;
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
                cout << "badly formatted nexus file, missing 'MATRIX' in data/character block" << endl;
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
        vector <string> tokens;
        string del(" \t");
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
        vector<string> tokens;
        string del(" \t");
        string tline;
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
            string tse = tokens[1];
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
        string curseq = "";
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
                    string id_ = tline.substr(1,tline.size()-1);
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
        string line1,line2,line3,line4;
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
vector <Sequence> read_interleaved_nexus_file (string filen, int ntax, int nchar) {
    vector <Sequence> seqs;
    string tline;
    ifstream infile(filen.c_str());
    bool done = false;
    
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
        cout << "badly formatted nexus file, missing 'MATRIX' in data/character block" << endl;
        exit(1);
    }
    
    int totcount = 0;
    int loopcount = 0;
    string del(" \t");
    while (getline(infile, tline)) {
        trim_spaces(tline);
        if (tline.size() != 0) {
            vector <string> tokens;
            tokenize(tline, tokens, del);
            if (tokens.size() > 1) {
                Sequence seq;
                for (unsigned int i=0; i < tokens.size(); i++) {
                    trim_spaces(tokens[i]);
                }
                if (tokens[0].compare(";") == 0) {
                    cout << "Huh?" << endl;
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
                string terp = seqs[ntax - 1].get_sequence();
                if (terp.size() == nchar) {
                    break;
                }
            }
        }
    }
    infile.close();
    //cout << "Seqs has " << seqs.size() << " taxa and "
    //        << seqs[0].get_sequence().size() << " sites." << endl;
    return seqs;
}

// don't search for MATRIX; if we know it is interleaved, MATRIX has already been read
vector <Sequence> read_interleaved_nexus (istream & stri, int ntax, int nchar) {
    vector <Sequence> seqs;
    string tline;
    bool done = false;
    
    int totcount = 0;
    int loopcount = 0;
    string del(" \t");
    while (getline(stri, tline)) {
        trim_spaces(tline);
        if (tline.size() != 0) {
            vector <string> tokens;
            tokenize(tline, tokens, del);
            if (tokens.size() > 1) {
                Sequence seq;
                for (unsigned int i=0; i < tokens.size(); i++) {
                    trim_spaces(tokens[i]);
                }
                if (tokens[0].compare(";") == 0) {
                    cout << "Huh?" << endl;
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
                string terp = seqs[ntax - 1].get_sequence();
                if (terp.size() == nchar) {
                    break;
                }
            }
        }
    }
    //cout << "Seqs has " << seqs.size() << " taxa and "
    //        << seqs[0].get_sequence().size() << " sites." << endl;
    return seqs;
}

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta)
 * TODO: need to add csv
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
 */
int test_char_filetype_stream(istream & stri, string & retstring) {
    if (!getline(stri, retstring)) {
        cout << "ERROR: end of file too soon" << endl;
    }
    int ret = 666; // if you get 666, there is no filetype set
    //NEXUS
    if (retstring[0] == '#') {
        ret = 0;
    } else if (retstring[0] == '>') {
        ret = 2;
    } else {    
        vector<string> tokens;
        string del(" \t");
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
bool read_next_seq_char_from_stream(istream & stri, int ftype, string & retstring, Sequence & seq) {
    string tline;
    if (ftype == 1) { // phylip
        vector<string> tokens;
        tokens.clear();
        string del(" \t");
        string tline;
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
        vector<string> tokens;
        string del(" \t");
        string tline;
        string curseq = "";
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
                    string id_ = tline.substr(1,tline.size()-1);
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

// file-version
// prolly get rid of this in favour of the stream-based one
void get_nexus_dimensions_file (string & filen, int & numTaxa, int & numChar, bool & interleave) {
    numTaxa = numChar = 0;
    string tline;
    string temp;
    ifstream infile(filen.c_str());
    while (getline(infile, tline)) {
        if (!tline.empty()) {
            // convert to uppercase
            tline = string_to_upper(tline);
            vector <string> searchtokens = tokenize(tline);
            if (searchtokens[0] == "DIMENSIONS") {
            // get rid of '=' and ';'. tokens then easy to deal with.
                replace(tline.begin(), tline.end(), '=', ' ');
                replace(tline.begin(), tline.end(), ';', ' ');
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
                replace(tline.begin(), tline.end(), '=', ' ');
                replace(tline.begin(), tline.end(), ';', ' ');
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
void get_nexus_dimensions (istream & stri, int & numTaxa, int & numChar, bool & interleave) {
    numTaxa = numChar = 0;
    string tline;
    //string temp;
    while (getline(stri, tline)) {
        if (!tline.empty()) {
            // convert to uppercase
            tline = string_to_upper(tline);
            vector <string> searchtokens = tokenize(tline);
            if (searchtokens[0] == "DIMENSIONS") {
            // get rid of '=' and ';'. tokens then easy to deal with.
                replace(tline.begin(), tline.end(), '=', ' ');
                replace(tline.begin(), tline.end(), ';', ' ');
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
                replace(tline.begin(), tline.end(), '=', ' ');
                replace(tline.begin(), tline.end(), ';', ' ');
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


// **************************** //
// *** Deprecated functions *** //
// **************************** //

// Leftover file-centric functions from before working with streams

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta), @ (fastq)
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
 */

int test_seq_filetype(string filen) {
    string tline;
    ifstream infile(filen.c_str());
    int ret = 666; // if you get 666, there is no filetype set
    while (getline(infile,tline)) {
        if (tline.size() < 1) {
            continue;
        }
        if (tline[0] == '#') {
            ret = 0;
            break;
        }
        vector<string> tokens;
        string del(" \t");
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
        cout << searchtokens[0] << " " << searchtokens[1] << endl;
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

