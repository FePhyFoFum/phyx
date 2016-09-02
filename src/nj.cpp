/*
 * NJOI.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: joe
 */

#include <vector>
#include <algorithm>
#include <map>

using namespace std;

#include "nj.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

/*Calculates the Q matrix
 * Original Matrix: All the Original Distances
 * ConvertedMatrix: Conversion to Q Matrix
 * LengthMatrix: Adjusted values for branch length calculations
 */
void NJOI::CalcQ(int const& NumbOfSequences, vector< vector<double> >& OriginalMatrix, 
    vector< vector<double> >& ConvertedMatrix, vector< vector<double> >& LengthMatrix) {

    ConvertedMatrix = OriginalMatrix;
    vector<double> Sums (NumbOfSequences, 0.0);
    
    for (int i = 0; i < NumbOfSequences; i++) {
        Sums[i] = sum(OriginalMatrix[i]);
        std::transform(ConvertedMatrix[i].begin(), ConvertedMatrix[i].end(), ConvertedMatrix[i].begin(), 
            std::bind1st(std::multiplies<double>(), (double) (NumbOfSequences - 2)));
    }
    
    for (int i = 0; i < NumbOfSequences; i++) {
        for (int j = 0; j < NumbOfSequences; j++) {
            if (i != j) {
                LengthMatrix[i][j] = abs(Sums[i] - Sums[j]);
                ConvertedMatrix[i][j] -= (Sums[i] + Sums[j]);
            }
        }
    }
}

/*Calculate the Branch Lengths
 * NewMatrix: Has original distances
 * LengthMatrix: Has adjusted lengths
 */
void NJOI::FetchLengths(int const& NumbOfSequences, vector< vector<double> > const& NewMatrix,
    vector< vector<double> >& LengthMatrix, int const& mini1, int const& mini2,
    double & brlength1, double & brlength2) {

    brlength1 = (NewMatrix[mini1][mini2] + (LengthMatrix[mini1][mini2] / (double)(NumbOfSequences - 2))) * 0.5;
    brlength2 = NewMatrix[mini1][mini2] - brlength1;

}

/*Tree_Update
 *Updates the Tree info and bypasses using a tree structure by storing parts in an array
 *NewMatrix: Has the adjusted values from the QMatrix calculation
*/
void NJOI::Tree_Update(string& newname, vector<string>& names, map<int, string>& NumbKeys,
    int& NumbOfSequences, vector< vector<double> >& NewMatrix, int& mini1, int& mini2,
    double& brlength1, double& brlength2) {
    
    int msize = (int)NewMatrix.size();
    
    //update the tree values, Tree Size is the node it is at
    vector<double> row_hits (NewMatrix[mini1]);
    vector<double> col_hits = (NewMatrix[mini2]);
    
    double ColRow = 0.0;
    double small_length = NewMatrix[mini1][mini2]; // neighbor based correction
    newname = "(" + names[mini1] + ":" + to_string(brlength2 / (double)nchar_) +  ","
        + names[mini2] + ":" + to_string(brlength1 / (double)nchar_) +  ")";
    
    // erase in backwards order as it preserves the indexes
    names.erase(names.begin()+mini2);
    names.erase(names.begin()+mini1);
    names.insert(names.begin(), newname);
    
    //Make Smaller Matrix
    vector< vector<double> > temp_matrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
    
    int count = 0;
    //Make a new First Row and Column
    for (int i = 0; i < msize; i++) {
        if (i != mini1 && i != mini2) {
            ColRow = (col_hits[i] + row_hits[i] - small_length) * 0.5;
            count++;
            temp_matrix[0][count] = ColRow;
            temp_matrix[count][0] = ColRow; // do we need top and bottom?
        }
    }
    
    //Need to fill the rest of the matrix up again
    int icount = 1;
    int jcount = 0;
    for (int i = 0; i < msize; i++) {
        jcount = 1;
        if (i != mini1 && i != mini2) {
            for (int j = 0; j < msize; j++) {
                if (j != mini1 && j != mini2) {
                    temp_matrix[icount][jcount] = NewMatrix[i][j];
                    jcount++;
                }
            }
            icount++;
        }
    }
    //cout << "NewMatrix.size() = " << NewMatrix.size() << "; temp_matrix.size() = "
    //    << temp_matrix.size() << "; NumbOfSequences = " << NumbOfSequences << endl;
    // NewMatrix is reduced b y1 in each dimension
    NewMatrix = temp_matrix;
}

// has to be a more efficient way of doing this!
void NJOI::Choose_Smallest(int& NumbOfSequences, vector< vector<double> > const& Matrix,
    int & mini1, int & mini2) {
    //super large value
    double MIN = 99999999999.99;
    for (int i = 0; i < (NumbOfSequences - 1); i++) {
        
        int idx = min_element(Matrix[i].begin() + (i + 1), Matrix[i].end()) - Matrix[i].begin();
        if (Matrix[i][idx] < MIN) {
            MIN = Matrix[i][idx];
            mini1 = i;
            mini2 = idx;
        }
    }
    NumbOfSequences--;
}

//NumbKeys Contains the names and their matching number
//Matrix Contains The original matrix
//Main Tree Making Matrix
void NJOI::TREEMAKE(vector<string>& names, map <int, string>& NumbKeys,
    vector< vector<double> >& Matrix) {

    int mini1 = 0, mini2 = 0;
    int NumbOfSequences = NumbKeys.size();
    double brlength1 = 0.0;
    double brlength2 = 0.0;
    vector< vector<double> > LengthMatrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
    string newname;
    while (NumbOfSequences > 2) {
        vector< vector<double> > QMatrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
        CalcQ(NumbOfSequences, Matrix, QMatrix, LengthMatrix);
        Choose_Smallest(NumbOfSequences, QMatrix, mini1, mini2);
        FetchLengths((NumbOfSequences + 1), Matrix, LengthMatrix, mini1, mini2,
            brlength1, brlength2);
        
        Tree_Update(newname, names, NumbKeys, NumbOfSequences, Matrix, mini1,
            mini2, brlength1, brlength2);
    }
    //double adjlength = (Matrix[mini1][mini2] / 2); // The final branch length
    double adjlength = (Matrix[mini1][mini2] / 2) / (double)nchar_;
    newname = "(" + names[mini1] + ":" + to_string(adjlength) +  "," + names[mini2]
        + ":" + to_string(adjlength) +  ")";
    newick_string_ = newname + ";";
}


vector< vector<double> > NJOI::BuildMatrix (map <string, string>& sequences) {

    vector<string> SequenceName;
    map <string, string>::iterator iter, iter2;
    string fasta, SeqName, MatchName;
    int FirstCount = 0;
    double MatchScore;

    // an easier way to initialize a vector of vectors:
    int ntax = sequences.size();
    vector< vector<double> > Score(ntax, vector<double>(ntax, 0.0));

    //compare all sequences to other sequences
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {
        fasta = iter -> second;
        SeqName = iter -> first;
        SequenceName.push_back(SeqName);
        int SecondCount = 0;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++) {
            MatchScore = (double) calc_hamming_dist(fasta, iter2 -> second);
            MatchName = SeqName + "," + iter2 -> first;
            Score[FirstCount][SecondCount] = MatchScore;
            SecondCount++;
        }
        FirstCount++;
    }
    return Score;
}



NJOI::NJOI (istream* pios, int & threads):ntax_(0), nchar_(0), nthreads_(threads) {
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    int seqcount = 0;
    // some error checking. should be in general seq reader class
    bool first = true;
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        sequences_[seq.get_id()] = seq.get_sequence();
        if (!first) {
            if ((int)seq.get_length() != nchar_) {
                cout << "Error: sequence " << seq.get_id() << " has "
                    << seq.get_length() << " characters, was expecting " 
                    << nchar_ << "." << endl << "Exiting." << endl;
                exit(1);
            }
        } else {
            nchar_ = seq.get_length();
            first = false;
        }
        seqcount++;
    }
    //fasta has a trailing one
    if (ft == 2) {
        sequences_[seq.get_id()] = seq.get_sequence();
        if ((int)seq.get_length() != nchar_) {
            cout << "Error: sequence " << seq.get_id() << " has "
                << seq.get_length() << " characters, was expecting " 
                << nchar_ << "." << endl << "Exiting." << endl;
            exit(1);
        };
        seqcount++;
    }
    ntax_ = seqcount;
    set_name_key ();
    Matrix = BuildMatrix(sequences_);
    TREEMAKE(names_, name_key_, Matrix);
}

void NJOI::set_name_key () {
    int count = 0;
    for (iter_ = sequences_.begin(); iter_ != sequences_.end(); iter_++) {
        name_key_[count] = iter_ -> first;
        names_.push_back(iter_ -> first);
        count++;
    }
}

string NJOI::get_newick () {
    return newick_string_;
}
