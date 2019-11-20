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
        sequences_[seq.get_id()] = seq.get_sequence();
        if (!first) {
            if ((int)seq.get_length() != nchar_) {
                std::cout << "Error: sequence " << seq.get_id() << " has "
                    << seq.get_length() << " characters, was expecting " 
                    << nchar_ << "." << std::endl << "Exiting." << std::endl;
                exit(1);
            }
        } else {
            nchar_ = seq.get_length();
            first = false;
        }
        nameKey_[seqcount] = seq.get_id();
        names_.push_back(seq.get_id());
        seqcount++;
    }
    //fasta has a trailing one
    if (ft == 2) {
        sequences_[seq.get_id()] = seq.get_sequence();
        if ((int)seq.get_length() != nchar_) {
            std::cout << "Error: sequence " << seq.get_id() << " has "
                << seq.get_length() << " characters, was expecting " 
                << nchar_ << "." << std::endl << "Exiting." << std::endl;
            exit(1);
        }
        nameKey_[seqcount] = seq.get_id();
        names_.push_back(seq.get_id());
        seqcount++;
    }
    ntax_ = seqcount;
    //set_name_key ();
    distmatrix_ = build_matrix(sequences_);
    make_tree(names_, nameKey_, distmatrix_);
}


void update_tree (std::string& newname, std::vector<std::string>& names, std::map<int, std::string>& numkeys, 
    int& node_list, std::vector< std::vector<double> >& newmatrix, int& mini1, int& mini2) {
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
    if (node_list > 1) {
        newname += "#" + std::to_string(br_length); // store height
    }

    names.erase(names.begin()+mini1);
    names.erase(names.begin()+(mini2-1));
    names.insert(names.begin(), newname);
    
    // Make Smaller Matrix
    std::vector< std::vector<double> > temp_matrix(node_list, std::vector<double>(node_list, 0.0));
    
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


void UPGMA::choose_small (int& node_list, const std::vector< std::vector<double> >& Matrix,
    int& mini1, int& mini2) {
    //super large value
    double MIN = 99999999999.99;
    for (int i = 0; i < (node_list - 1); i++) {
        int idx = std::min_element(Matrix[i].begin() + (i + 1), Matrix[i].end()) - Matrix[i].begin();
        if (Matrix[i][idx] < MIN) {
            MIN = Matrix[i][idx];
            mini1 = i;
            mini2 = idx;
        }
    }
    node_list--;
}


// numkeys Contains the names and their matching number
// Matrix Contains The original matrix
void UPGMA::make_tree (std::vector<std::string>& names, std::map<int, std::string>& numkeys,
    std::vector< std::vector<double> >& Matrix) {

    int mini1 = 0, mini2 = 0;
    int NumbOfSequences = numkeys.size();
    std::vector< std::vector<double> > newmatrix;
    std::map<int, std::string>::iterator iter;
    std::string newname;
    while (NumbOfSequences > 1) {
        choose_small(NumbOfSequences, Matrix, mini1, mini2);
        update_tree(newname, names, numkeys, NumbOfSequences, Matrix, mini1, mini2);
    }
    newickstring_ = newname + ";";
}


std::vector< std::vector<double> > UPGMA::build_matrix (std::map<std::string, std::string>& sequences) {

    std::vector<std::string> SequenceName;
    std::map<std::string, std::string>::iterator iter, iter2;
    std::string fasta, SeqName, MatchName;
    int FirstCount = 0;
    double MatchScore;

    // an easier way to initialize a vector of vectors:
    std::vector< std::vector<double> > Score(ntax_, std::vector<double>(ntax_, 0.0));
    
    //compare all sequences to other sequences
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {
        fasta = iter -> second;
        SeqName = iter -> first;
        SequenceName.push_back(SeqName);
        int SecondCount = 0;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++) {
            //MatchScore = CalcSeqDiffs(fasta, iter2 -> second);
            MatchScore = (double) calc_hamming_dist(fasta, iter2 -> second);
            MatchName = SeqName + "," + iter2 -> first;
            Score[FirstCount][SecondCount] = MatchScore;
            SecondCount++;
        }
        FirstCount++;

    }
    //prints the distance matrix maybe too verbose
    std::cout << "\t";
    for (unsigned int i = 0; i < SequenceName.size(); i++) {
        std::cout << SequenceName[i] << "\t";
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < Score.size(); i++) {
        std::cout << SequenceName[i] << "\t";
        for (unsigned int j = 0; j < Score[i].size(); j++) {
            std::cout << Score[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    return Score;
}


// populate these when reading in the sequences
void UPGMA::set_name_key () {
    int count = 0;
    std::map<std::string, std::string>::iterator iter;
    for (iter = sequences_.begin(); iter != sequences_.end(); iter++) {
        nameKey_[count] = iter -> first;
        names_.push_back(iter -> first);
        count++;
    }
}


std::string UPGMA::get_newick () const {
    return newickstring_;
}
