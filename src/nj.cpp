/*
 * NJOI.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: joe
 */

#include <string>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
#include <map>
#include <iterator>

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
void CalcQ(int const& NumbOfSequences, vector< vector<double> >& OriginalMatrix, 
    vector< vector<double> >& ConvertedMatrix, vector< vector<double> >& LengthMatrix) {

    ConvertedMatrix = OriginalMatrix;
//    double row_sum = 0.0;
//    double col_sum = 0.0;
//    vector<double> Sums ((NumbOfSequences + 1), 0.0);
    vector<double> Sums (NumbOfSequences, 0.0);
    
    //double xx = (double) NumbOfSequences - 2.0;
    
    for (int i = 0; i < NumbOfSequences; i++) {
//        row_sum = 0.0;
//        for (int j = 0; j < NumbOfSequences; j++) {
//            row_sum += ConvertedMatrix[i][j];
//            ConvertedMatrix[i][j] *= (NumbOfSequences - 2);
//        }
//        if (i == 0) {
//            col_sum = row_sum;
//        }
        //Sums[i] = calculate_vector_double_sum(OriginalMatrix[i]);
        Sums[i] = sum(OriginalMatrix[i]);
        std::transform(ConvertedMatrix[i].begin(), ConvertedMatrix[i].end(), ConvertedMatrix[i].begin(), 
            std::bind1st(std::multiplies<double>(), (double) (NumbOfSequences - 2)));
    }
//    col_sum = Sums[0];
//    Sums[NumbOfSequences] = Sums[0];
    //Sums.push_back(col_sum); // so first and last value of Sums are the same? wait, never gets this big...
    
    for (int i = 0; i < NumbOfSequences; i++) {
        for (int j = 0; j < NumbOfSequences; j++) {
            if (i != j) {
                LengthMatrix[i][j] = abs(Sums[i] - Sums[j]);
                ConvertedMatrix[i][j] -= (Sums[i] + Sums[j]);
//                ConvertedMatrix[i][j] -= Sums[i];
//                ConvertedMatrix[i][j] -= Sums[j];
            }
        }
    }
    //Sums.clear();
}
/*Calculate the Branch Lengths
 * NewMatrix: Has original distances
 * LengthMatrix: Has adjusted lengths
 */
void FetchLengths(int const& NumbOfSequences, vector< vector<double> > const& NewMatrix,
    vector< vector<double> >& LengthMatrix, int const& mini1, int const& mini2,
    double & brlength1, double & brlength2) {

    //Cant figure out how to bypass changing to a string then a double, but their has to
    //be a command to change integer to double
    //string hold = to_string(NumbOfSequences);
    //brlength1 = (NewMatrix[mini1][mini2] / 2.0) + (1 / (2*(stod(hold)-2))) * LengthMatrix[mini1][mini2];
    brlength1 = (NewMatrix[mini1][mini2] + (LengthMatrix[mini1][mini2] / (double)(NumbOfSequences - 2))) * 0.5;
    brlength2 = NewMatrix[mini1][mini2] - brlength1;

}

/*Tree_Update
 *Updates the Tree info and bypasses using a tree structure by storing parts in an array
 *NewMatrix: Has the adjusted values from the QMatrix calculation
 */
// NumbOfSequences arg wasn't being used
//void Tree_Update(int NumbOfSequences, string& newname, vector<string>& names,
//    map<int, string>& NumbKeys, int& node_list, vector< vector<double> >& NewMatrix,
//    int& mini1, int& mini2, double& brlength1, double& brlength2) {
void Tree_Update(string& newname, vector<string>& names, map<int, string>& NumbKeys,
    int& node_list, vector< vector<double> >& NewMatrix, int& mini1, int& mini2,
    double& brlength1, double& brlength2) {
    
    int msize = (int)NewMatrix.size();
    
    //update the tree values, Tree Size is the node it is at
    //vector<double> row_hits (msize, 0.0); // this is simply NewMatrix[mini1]
    vector<double> row_hits (NewMatrix[mini1]);
    //vector<double> col_hits (msize, 0.0); // this is simply NewMatrix[mini1]
    vector<double> col_hits = (NewMatrix[mini2]);
    //vector<double> new_ColRow (msize, 0.0);
    //vector<double> row_hits, col_hits, new_ColRow;
    
    double ColRow = 0.0;
    double small_length = NewMatrix[mini1][mini2]; // neighbor based correction
    newname = "(" + names[mini1] + ":" + to_string(brlength2) +  ","
        + names[mini2] + ":" + to_string(brlength1) +  ")";
    
    // erase in backwards order as it preserves the indexes
    names.erase(names.begin()+mini2);
    names.erase(names.begin()+mini1);
    //names.erase(names.begin()+(mini2-1));
    names.insert(names.begin(), newname);
    
    //Make Smaller Matrix
    vector< vector<double> > temp_matrix(node_list, vector<double>(node_list, 0.0));
    
    //Reformat Matrix Takes in the hits
/*
    int rcount = 0, ccount = 0;
    for (int i = 0; i < msize; i++) {
        for (int j = 0; j < msize; j++) {
            if (i == mini1) {
                //row_hits.push_back(NewMatrix[mini1][j]);
                row_hits[rcount] = NewMatrix[mini1][j];
                rcount++;
            } else
            if (i == mini2) {
                //col_hits.push_back(NewMatrix[mini2][j]);
                col_hits[ccount] = NewMatrix[mini2][j];
                ccount++;
            }
        }
    }
*/
    
    int count = 0;
    //Make a new First Row and Column
    for (int i = 0; i < msize; i++) {
        //ColRow = (col_hits[i] + row_hits[i] - small_length) / 2;
        if (i != mini1 && i != mini2) {
            ColRow = (col_hits[i] + row_hits[i] - small_length) * 0.5;
//            if (i != mini2) {
                count++;
                temp_matrix[0][count] = ColRow;
                temp_matrix[count][0] = ColRow; // do we need top and bottom?
//            }
        }
        //new_ColRow.push_back(ColRow);
        //new_ColRow[count] = ColRow; // what is this for? not used
    }
    
    //Need to fill the rest of the matrix up again
    int icount = 1;
    int jcount = 0;
    for (int i = 0; i < msize; i++) {
        jcount = 1;
        if (i != mini1 && i != mini2) {
//            if(i != mini2) {
            for (int j = 0; j < msize; j++) {
                if (j != mini1 && j != mini2) {
//                    if (j != mini2) {
                        temp_matrix[icount][jcount] = NewMatrix[i][j];
                        jcount++;
//                    }
                }
            }
            icount++;
//            }
        }
    }
    //cout << "NewMatrix.size() = " << NewMatrix.size() << "; temp_matrix.size() = "
    //    << temp_matrix.size() << "; node_list = " << node_list << endl;
    // NewMatrix is reduced b y1 in each dimension
    NewMatrix = temp_matrix;
}

// has to be a more efficient way of doing this!
void Choose_Smallest(int& node_list, vector< vector<double> > const& Matrix,
    int & mini1, int & mini2) {
    //super large value
    double MIN = 99999999999.99;
    for (int i = 0; i < (node_list - 1); i++) {
        
        int idx = min_element(Matrix[i].begin() + (i + 1), Matrix[i].end()) - Matrix[i].begin();
        if (Matrix[i][idx] < MIN) {
            MIN = Matrix[i][idx];
            mini1 = i;
            mini2 = idx;
        }
        
//        for(int j = i + 1; j < node_list; j++) {
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
//Main Tree Making Matrix
void NJOI::TREEMAKE(vector<string>& names, map <int, string>& NumbKeys,
    vector< vector<double> >& Matrix) {

    int mini1 = 0, mini2 = 0;
    int NumbOfSequences = NumbKeys.size();
    double brlength1 = 0.0;
    double brlength2 = 0.0;
//    vector< vector<double> > NewMatrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
    vector< vector<double> > LengthMatrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
//    NewMatrix = Matrix; // not used
//    map<int, string>::iterator iter; // not used
    string newname;
    while (NumbOfSequences > 2) {
        vector< vector<double> > QMatrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
        CalcQ(NumbOfSequences, Matrix, QMatrix, LengthMatrix);
        Choose_Smallest(NumbOfSequences, QMatrix, mini1, mini2);
        FetchLengths((NumbOfSequences + 1), Matrix, LengthMatrix, mini1, mini2,
            brlength1, brlength2);
        
        // first arg wasn't being used
//        Tree_Update((NumbOfSequences + 1), newname, names, NumbKeys,
//            NumbOfSequences, Matrix, mini1, mini2, brlength1, brlength2);
        Tree_Update(newname, names, NumbKeys, NumbOfSequences, Matrix, mini1,
            mini2, brlength1, brlength2);
    }
    double adjlength = (Matrix[mini1][mini2] / 2); // The final branch length
    //cout << Matrix[mini1][mini2] << endl;
    //FetchLengths((NumbOfSequences + 1), Matrix, LengthMatrix, mini1, mini2, brlength1, brlength2);
    newname = "(" + names[mini1] + ":" + to_string(adjlength) +  "," + names[mini2] + ":" + to_string(adjlength) +  ")";
    newickstring = newname + ";";
}

// not used. doesn't have a declaration either
//double CalcSeqsDiffs(string const& sequence1, string const& sequence2) {
//    double score = 0;
//    for (unsigned int i = 0; i < sequence1.size(); i++) {
//        if (sequence1[i] != sequence2[i]) {
//            score++;
//        }
//    }
//    return score;
//}

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
            //MatchScore = CalcSeqsDiffs(fasta, iter2 -> second);
            MatchScore = (double) calc_hamming_dist(fasta, iter2 -> second);
            MatchName = SeqName + "," + iter2 -> first;
            Score[FirstCount][SecondCount] = MatchScore;
            SecondCount++;
        }
        FirstCount++;
    }
    return Score;
}


NJOI::NJOI() {
    // TODO Auto-generated constructor stub

}

NJOI::NJOI (istream* pios, int & threads):ntax(0), nchar(0), nthreads(threads) {
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    int seqcount = 0;
    // some error checking. should be in general seq reader class
    bool first = true;
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        sequences[seq.get_id()] = seq.get_sequence();
        if (!first) {
            if ((int)seq.get_length() != nchar) {
                cout << "Error: sequence " << seq.get_id() << " has "
                    << seq.get_length() << " characters, was expecting " 
                    << nchar << "." << endl << "Exiting." << endl;
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
            cout << "Error: sequence " << seq.get_id() << " has "
                << seq.get_length() << " characters, was expecting " 
                << nchar << "." << endl << "Exiting." << endl;
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
void NJOI::set_name_key () {
    int count = 0;
    for(iter = sequences.begin(); iter != sequences.end(); iter++) {
        NameKey[count] = iter -> first;
        names.push_back(iter -> first);
        count++;
    }
}

string NJOI::get_newick () {
    return newickstring;
}
