#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include "upgma.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"


UPGMA::UPGMA (std::istream* pios):ntax_(0), nchar_(0) {
    Sequence seq;
    std::string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    int seqcount = 0;
    // some error checking. should be in general seq reader class
    bool first = true;
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        seqs_.push_back(seq);
        
        sequences_[seq.get_id()] = seq.get_sequence();
        if (!first) {
            if ((int)seq.get_length() != nchar_) {
                std::cerr << "Error: sequence " << seq.get_id() << " has "
                    << seq.get_length() << " characters, was expecting " 
                    << nchar_ << ". Exiting." << std::endl;
                exit(1);
            }
        } else {
            nchar_ = seq.get_length();
            first = false;
        }
        //nameKey_[seqcount] = seq.get_id();
        names_.push_back(seq.get_id());
        seqcount++;
    }
    //seq has a trailing one
    if (ft == 2) {
        seqs_.push_back(seq);
        
        sequences_[seq.get_id()] = seq.get_sequence();
        if ((int)seq.get_length() != nchar_) {
            std::cerr << "Error: sequence " << seq.get_id() << " has "
                << seq.get_length() << " characters, was expecting " 
                << nchar_ << ". Exiting." << std::endl;
            exit(1);
        }
        //nameKey_[seqcount] = seq.get_id();
        names_.push_back(seq.get_id());
        seqcount++;
    }
    ntax_ = seqcount;
    //distmatrix_ = build_matrix(sequences_);
    distmatrix_ = build_matrix();
    //make_tree(names_, nameKey_, distmatrix_);
    make_tree();
}


// get rid of all the map stuff
//std::vector< std::vector<double> > UPGMA::build_matrix (std::map<std::string, std::string>& sequences) {
std::vector< std::vector<double> > UPGMA::build_matrix () {

    /*
    //std::vector<std::string> SequenceName; // but this is already present as names_
    std::map<std::string, std::string>::iterator iter, iter2;
    std::string seq, SeqName, MatchName;
    int FirstCount = 0;
    double MatchScore;

    // an easier way to initialize a vector of vectors:
    std::vector< std::vector<double> > Score(ntax_, std::vector<double>(ntax_, 0.0));
    
    // compare all sequences to other sequences
    // really only need half the matrix
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {
        seq = iter -> second;
        SeqName = iter -> first;
        //SequenceName.push_back(SeqName);
        int SecondCount = 0;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++) {
            //MatchScore = CalcSeqDiffs(seq, iter2 -> second);
            MatchScore = (double) calc_hamming_dist(seq, iter2 -> second);
            MatchName = SeqName + "," + iter2 -> first;
            //std::cout << "MatchName = " << MatchName << std::endl;
            Score[FirstCount][SecondCount] = MatchScore;
            SecondCount++;
        }
        FirstCount++;
    }
    */
    
    // smarter version:
    // 1) skip self comparisons
    // 2) only calculate one half of matrix (i.e., no duplicate calcs)
    std::vector< std::vector<double> > distances(ntax_, std::vector<double>(ntax_, 0.0));
    double tempScore = 0.0;
    for (int i = 0; i < ntax_; i++) {
        std::string seq1 = seqs_[i].get_sequence();
        for (int j = (i + 1); j < ntax_; j++) {
            std::string seq2 = seqs_[j].get_sequence();
            // get distance
            tempScore = (double) calc_hamming_dist(seq1, seq2);
            // put in both top and bottom of matrix
            distances[i][j] = distances[j][i] = tempScore;
        }
    }
    
    // don't want this
    /*
    std::cout << std::endl << "OLD" << std::endl;
    std::cout << "\t";
    for (unsigned int i = 0; i < names_.size(); i++) {
        std::cout << names_[i] << "\t";
    }
    std::cout << std::endl;
    for (int i = 0; i < ntax_; i++) {
        std::cout << names_[i] << "\t";
        for (int j = 0; j < ntax_; j++) {
            std::cout << Score[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    */
    //std::cout << std::endl << "NEW" << std::endl;
    
    std::cout << "\t";
    for (unsigned int i = 0; i < names_.size(); i++) {
        std::cout << names_[i] << "\t";
    }
    std::cout << std::endl;
    for (int i = 0; i < ntax_; i++) {
        std::cout << names_[i] << "\t";
        for (int j = 0; j < ntax_; j++) {
            std::cout << distances[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    
    return distances;
}


void update_tree (std::string& newname, std::vector<std::string>& names,
    int& numseqs, std::vector< std::vector<double> >& newmatrix, int& mini1, int& mini2) {
    //update the tree values, Tree Size is the node it is at
    std::vector<double> row_hits, col_hits, new_ColRow;
    double br_length = newmatrix[mini1][mini2] / 2.0;
        
    std::string length1 = std::to_string(br_length);
    std::string length2 = length1;
    double ColRow = 0.0;
    int matrixsize = newmatrix.size();

    // extremely hacky way to get correct ELs
    std::size_t found = names[mini1].find("#");
    if (found != std::string::npos) {
        std::string terp = names[mini1];
        std::size_t pos = terp.find("#"); 
        double oldheight = std::stod(terp.substr(pos+1));
        names[mini1] = terp.substr(0, pos);
        length1 = std::to_string(br_length - oldheight);
    }
    // have to do it for the other side too
    found = names[mini2].find("#");
    if (found != std::string::npos) {
        std::string terp = names[mini2];
        std::size_t pos = terp.find("#"); 
        double oldheight = std::stod(terp.substr(pos+1));
        names[mini2] = terp.substr(0, pos);
        length2 = std::to_string(br_length - oldheight);
    }

    newname = "(" + names[mini1] + ":" + length1 + "," + names[mini2] + ":" + length2 + ")";
    if (numseqs > 1) {
        newname += "#" + std::to_string(br_length); // store height
    }

    names.erase(names.begin()+mini1);
    names.erase(names.begin()+(mini2-1));
    names.insert(names.begin(), newname);
    
    // Make Smaller Matrix
    std::vector< std::vector<double> > temp_matrix(numseqs, std::vector<double>(numseqs, 0.0));
    
    // Reformat Matrix
    for (int i = 0; i < matrixsize; i++) {
        for (int j = 0; j < matrixsize; j++) {
            if (i == mini1) {
                row_hits.push_back(newmatrix[mini1][j]);
            } else if (i == mini2) {
                col_hits.push_back(newmatrix[mini2][j]);
            }
        }
    }
    
    int count = 0;
    //Make a new First Row and Column
    for (int i = 0; i < (int)col_hits.size(); i++) {
        ColRow = (col_hits[i] + row_hits[i]) / 2;
        if (i != mini1) {
            if (i != mini2) {
                count++;
                temp_matrix[0][count] = ColRow;
                temp_matrix[count][0] = ColRow;
            }
        }
        new_ColRow.push_back(ColRow);
    }
    
    // Need to fill the rest of the matrix up again
    int icount = 1;
    int jcount = 0;
    //std::cout << "newmatrix.size() = " << newmatrix.size()
    //    << "; mini1 = " << mini1 << "; mini2 = " << mini2 << std::endl;
    // mini1 is always < mini2
    for (int i = 0; i < matrixsize; i++) {
        //print_vector(newmatrix[i]);
        //std::copy(newmatrix[i].begin(), newmatrix[i].end(), std::ostream_iterator<double>(std::cout, " "));
        jcount = 1;
        if (i != mini1 && i != mini2) {
//            if (i != mini2) {
            for (int j = 0; j < matrixsize; j++) {
                if (j != mini1 && j != mini2) {
                    //if (j != mini2) {
                        temp_matrix[icount][jcount] = newmatrix[i][j];
                        jcount++;
                    //}
                }
            }
            icount++;
//            }
        }
    }
    newmatrix = temp_matrix;
}


// find smallest pairwise distance. will always find this on the top half of the matrix
void UPGMA::choose_small (int& numseqs, const std::vector< std::vector<double> >& dmatrix,
    int& mini1, int& mini2) {
    //super large value
    double MIN = 99999999999.99;
    for (int i = 0; i < (numseqs - 1); i++) {
        int idx = std::min_element(dmatrix[i].begin() + (i + 1), dmatrix[i].end()) - dmatrix[i].begin();
        if (dmatrix[i][idx] < MIN) {
            MIN = dmatrix[i][idx];
            mini1 = i;
            mini2 = idx;
        }
    }
    numseqs--;
}


// numkeys contains the names and their matching number. this is never used
// Matrix contains the original matrix
void UPGMA::make_tree () {
    // initialize
    int mini1 = 0, mini2 = 0;
    std::vector<std::string> names = names_;
    std::vector< std::vector<double> > dmatrix = distmatrix_;
    int numseqs = ntax_;
    std::string newname;
    
    while (numseqs > 1) {
        choose_small(numseqs, dmatrix, mini1, mini2);
        //update_tree(newname, names, numkeys, numseqs, Matrix, mini1, mini2);
        update_tree(newname, names, numseqs, dmatrix, mini1, mini2);
    }
    newickstring_ = newname + ";";
}


std::string UPGMA::get_newick () const {
    return newickstring_;
}
