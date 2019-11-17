#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <map>
#include <iterator>
#include <cstring>
#include <getopt.h>

#include "upgma.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

void update_tree(std::string& newname, std::vector<std::string>& names, std::map<int, std::string>& NumbKeys, 
    int& node_list, std::vector< std::vector<double> >& NewMatrix, int& mini1, int& mini2) {
    //update the tree values, Tree Size is the node it is at
    std::vector<double> row_hits, col_hits, new_ColRow;
    double br_length = NewMatrix[mini1][mini2] / 2.0;
        
    std::string length1 = std::to_string(br_length);
    std::string length2 = length1;
    double ColRow = 0.0;
    int matrixsize = NewMatrix.size();

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
    
    //Make Smaller Matrix
    std::vector< std::vector<double> > temp_matrix(node_list, std::vector<double>(node_list, 0.0));
//    vector< vector<double> > temp_matrix;    
//    for (int i = 0; i < node_list; i++) {
//        vector<double> row;
//        for (int j = 0; j < node_list; j++) {
//            row.push_back(0);
//        }
//        temp_matrix.push_back(row);
//    }
    
    // Reformat Matrix
    for (int i = 0; i < matrixsize; i++) {
        for (int j = 0; j < matrixsize; j++) {
            if (i == mini1) {
                row_hits.push_back(NewMatrix[mini1][j]);
            } else if (i == mini2) {
                col_hits.push_back(NewMatrix[mini2][j]);
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
    
    //Need to fill the rest of the matrix up again
    int icount = 1;
    int jcount = 0;
    // std::cout << "NewMatrix.size() = " << NewMatrix.size()
    //    << "; mini1 = " << mini1 << "; mini2 = " << mini2 << std::endl;
    // mini1 is always < mini2
    for (int i = 0; i < matrixsize; i++) {
        // print vectorcontents
        //print_vector(NewMatrix[i]);
        //std::copy(NewMatrix[i].begin(), NewMatrix[i].end(), std::ostream_iterator<double>(std:: std::cout, " "));
        jcount = 1;
        if (i != mini1 && i != mini2) {
//            if (i != mini2) {
            for (int j = 0; j < matrixsize; j++) {
                if (j != mini1 && j != mini2) {
                    //if (j != mini2) {
                        temp_matrix[icount][jcount] = NewMatrix[i][j];
                        jcount++;
                    //}
                }
            }
            icount++;
//            }
        }
    }
    //Print the new matrix makes it very verbose but
    //fun to watch
    /*
     std::cout << name1 << "  <===============NewMatrix==============>  " << name2 << std::endl;
    for (int i = 0; i < temp_matrix.size(); i++) {
        for (int j = 0; j < temp_matrix.size(); j++) {
             std::cout <<temp_matrix[i][j] << "\t";
        }
         std::cout << std::endl;
    }*/
    NewMatrix = temp_matrix;
}


void Choose_Small(int& node_list, std::vector< std::vector<double> > const& Matrix,
    int& mini1, int& mini2) {
    //super large value
    double MIN = 99999999999.99;
    for (int i = 0; i < (node_list - 1); i++) {
        
        int idx = min_element(Matrix[i].begin() + (i + 1), Matrix[i].end()) - Matrix[i].begin();
        if (Matrix[i][idx] < MIN) {
            MIN = Matrix[i][idx];
            mini1 = i;
            mini2 = idx;
        }
        
//        for (int j = i + 1; j < node_list; j++) {
//            if (Matrix[i][j] < MIN) {
//                MIN = Matrix[i][j];
//                mini1 = i;
//                mini2 = j;
//            }
//        }
    }
    node_list--;
}


//NumbKeys Contains the names and their matching number
//Matrix Contains The original matrix
void UPGMA::TREEMAKE(std::vector<std::string>& names, std::map<int, std::string>& NumbKeys,
    std::vector< std::vector<double> >& Matrix) {

    int mini1 = 0, mini2 = 0;
    int NumbOfSequences = NumbKeys.size();
    std::vector< std::vector<double> > NewMatrix;
    std::map<int, std::string>::iterator iter;
    std::string newname;
    while (NumbOfSequences > 1) {
        Choose_Small(NumbOfSequences, Matrix, mini1, mini2);
        update_tree(newname, names, NumbKeys, NumbOfSequences, Matrix, mini1, mini2);
    }
    newickstring = newname + ";";
}


// not used
double CalcSeqDiffs(std::string& sequence1, std::string& sequence2) {
    double score = 0;
    for (unsigned int i = 0; i < sequence1.size(); i++) {
        if (sequence1[i] != sequence2[i]) {
            score++;
        }
    }
    return score;
}


std::vector< std::vector<double> > UPGMA::BuildMatrix (std::map<std::string, std::string>& sequences) {

    std::vector<std::string> SequenceName;
    std::map<std::string, std::string>::iterator iter, iter2;
    std::string fasta, SeqName, MatchName;
    int FirstCount = 0;
    double MatchScore;

    // an easier way to initialize a vectorof vectors:
    int ntax = sequences.size(); // this should be a property of the class, calculated once
    std::vector< std::vector<double> > Score(ntax, std::vector<double>(ntax, 0.0));
    
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


UPGMA::UPGMA() {
    // TODO Auto-generated constructor stub
}


// *** some alternate functions below *** //


// alternate constructor
UPGMA::UPGMA (std::istream* pios):ntax(0), nchar(0) {
    Sequence seq;
    std::string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    int seqcount = 0;
    // some error checking. should be in general seq reader class
    bool first = true;
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        sequences[seq.get_id()] = seq.get_sequence();
        if (!first) {
            if ((int)seq.get_length() != nchar) {
                 std::cout << "Error: sequence " << seq.get_id() << " has "
                    << seq.get_length() << " characters, was expecting " 
                    << nchar << "." << std::endl << "Exiting." << std::endl;
                exit(1);
            }
        } else {
            nchar = seq.get_length();
            first = false;
        }
        NameKey[seqcount] = seq.get_id();
        names.push_back(seq.get_id());
        seqcount++;
    }
    //fasta has a trailing one
    if (ft == 2) {
        sequences[seq.get_id()] = seq.get_sequence();
        if ((int)seq.get_length() != nchar) {
             std::cout << "Error: sequence " << seq.get_id() << " has "
                << seq.get_length() << " characters, was expecting " 
                << nchar << "." << std::endl << "Exiting." << std::endl;
            exit(1);
        }
        NameKey[seqcount] = seq.get_id();
        names.push_back(seq.get_id());
        seqcount++;
    }
    ntax = seqcount;
    //set_name_key ();
    Matrix = BuildMatrix(sequences);
    TREEMAKE(names, NameKey, Matrix);
}


// populate these when reading in the sequences
void UPGMA::set_name_key () {
    int count = 0;
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {
        NameKey[count] = iter -> first;
        names.push_back(iter -> first);
        count++;
    }
}


std::string UPGMA::get_newick () {
    return newickstring;
}
