#include <vector>
#include <algorithm>
#include <map>

using namespace std; // commenting this out changes the result?!? ***

#include "nj.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"


/*Calculates the Q matrix
 * Original Matrix: All the Original Distances
 * ConvertedMatrix: Conversion to Q Matrix
 * LengthMatrix: Adjusted values for branch length calculations
 */
void NJOI::CalcQ (const int& NumbOfSequences, std::vector< std::vector<double> >& OriginalMatrix, 
    std::vector< std::vector<double> >& ConvertedMatrix, std::vector< std::vector<double> >& LengthMatrix) {

    ConvertedMatrix = OriginalMatrix;
    std::vector<double> Sums(NumbOfSequences, 0.0);
    
    for (int i = 0; i < NumbOfSequences; i++) {
        Sums[i] = sum(OriginalMatrix[i]);
        std::transform(ConvertedMatrix[i].begin(), ConvertedMatrix[i].end(), ConvertedMatrix[i].begin(), 
            std::bind1st(std::multiplies<double>(), static_cast<double>(NumbOfSequences - 2)));
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
void NJOI::FetchLengths (const int& NumbOfSequences, const std::vector< std::vector<double> >& NewMatrix,
    std::vector< std::vector<double> >& LengthMatrix, const int& mini1, const int& mini2,
    double & brlength1, double & brlength2) {

    brlength1 = (NewMatrix[mini1][mini2] + (LengthMatrix[mini1][mini2] / static_cast<double>(NumbOfSequences - 2))) * 0.5;
    brlength2 = NewMatrix[mini1][mini2] - brlength1;
}


/*Tree_Update
 *Updates the Tree info and bypasses using a tree structure by storing parts in an array
 *NewMatrix: Has the adjusted values from the QMatrix calculation
*/
void NJOI::Tree_Update (std::string& newname, std::vector<std::string>& names,
    int& NumbOfSequences, std::vector< std::vector<double> >& NewMatrix, int& mini1, int& mini2,
    double& brlength1, double& brlength2) {
    
    int msize = static_cast<int>(NewMatrix.size());
    
    //update the tree values, Tree Size is the node it is at
    std::vector<double> row_hits = NewMatrix[mini1];
    std::vector<double> col_hits = NewMatrix[mini2];
    
    double ColRow = 0.0;
    double small_length = NewMatrix[mini1][mini2]; // neighbor based correction
    newname = "(" + names[mini1] + ":" + std::to_string(brlength2 / static_cast<double>(num_char_)) +  ","
        + names[mini2] + ":" + std::to_string(brlength1 / static_cast<double>(num_char_)) +  ")";
    
    // erase in backwards order as it preserves the indexes
    names.erase(names.begin()+mini2);
    names.erase(names.begin()+mini1);
    names.insert(names.begin(), newname);
    
    // Make Smaller Matrix
    std::vector< std::vector<double> > temp_matrix(NumbOfSequences,
            std::vector<double>(NumbOfSequences, 0.0));
    
    int count = 0;
    // Make a new First Row and Column
    for (int i = 0; i < msize; i++) {
        if (i != mini1 && i != mini2) {
            ColRow = (col_hits[i] + row_hits[i] - small_length) * 0.5;
            count++;
            temp_matrix[0][count] = ColRow;
            temp_matrix[count][0] = ColRow; // do we need top and bottom?
        }
    }
    
    // Need to fill the rest of the matrix up again
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
    //std::cout << "NewMatrix.size() = " << NewMatrix.size() << "; temp_matrix.size() = "
    //    << temp_matrix.size() << "; NumbOfSequences = " << NumbOfSequences << std::endl;
    // NewMatrix is reduced by 1 in each dimension
    NewMatrix = temp_matrix;
}


// has to be a more efficient way of doing this!
void NJOI::Choose_Smallest (int& NumbOfSequences, const std::vector< std::vector<double> >& Matrix,
    int & mini1, int & mini2) {
    //super large value
    double MIN = 99999999999.99;
    for (int i = 0; i < (NumbOfSequences - 1); i++) {
        
        int idx = std::min_element(Matrix[i].begin() + (i + 1), Matrix[i].end()) - Matrix[i].begin();
        if (Matrix[i][idx] < MIN) {
            MIN = Matrix[i][idx];
            mini1 = i;
            mini2 = idx;
        }
    }
    NumbOfSequences--;
}


// NumbKeys Contains the names and their matching number
// Matrix Contains The original matrix
// Main Tree Making Matrix
void NJOI::TREEMAKE (std::vector<std::string>& names, std::map<int, std::string>& NumbKeys,
    std::vector< std::vector<double> >& Matrix) {

    int mini1 = 0, mini2 = 0;
    int NumbOfSequences = NumbKeys.size();
    double brlength1 = 0.0;
    double brlength2 = 0.0;
    std::vector< std::vector<double> > LengthMatrix(NumbOfSequences, std::vector<double>(NumbOfSequences, 0.0));
    std::string newname;
    while (NumbOfSequences > 2) {
        std::vector< std::vector<double> > QMatrix(NumbOfSequences, std::vector<double>(NumbOfSequences, 0.0));
        CalcQ(NumbOfSequences, Matrix, QMatrix, LengthMatrix);
        Choose_Smallest(NumbOfSequences, QMatrix, mini1, mini2);
        FetchLengths((NumbOfSequences + 1), Matrix, LengthMatrix, mini1, mini2,
            brlength1, brlength2);
        
        Tree_Update(newname, names, NumbOfSequences, Matrix, mini1,
            mini2, brlength1, brlength2);
    }
    //double adjlength = (Matrix[mini1][mini2] / 2); // The final branch length
    double adjlength = (Matrix[mini1][mini2] / 2) / static_cast<double>(num_char_);
    newname = "(" + names[mini1] + ":" + std::to_string(adjlength) +  "," + names[mini2]
        + ":" + std::to_string(adjlength) +  ")";
    newick_string_ = newname + ";";
}


std::vector< std::vector<double> > NJOI::BuildMatrix (std::map<std::string, std::string>& sequences) {

    std::vector<std::string> SequenceName;
    std::map <std::string, std::string>::iterator iter, iter2;
    std::string fasta, SeqName, MatchName;
    int FirstCount = 0;
    double MatchScore;

    // an easier way to initialize a std::vector of std::vectors:
    int ntax = sequences.size();
    std::vector< std::vector<double> > Score(ntax, std::vector<double>(ntax, 0.0));

    //compare all sequences to other sequences
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {
        fasta = iter -> second;
        SeqName = iter -> first;
        SequenceName.push_back(SeqName);
        int SecondCount = 0;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++) {
            MatchScore = static_cast<double>(calc_hamming_dist(fasta, iter2 -> second));
            MatchName = SeqName + "," + iter2 -> first;
            Score[FirstCount][SecondCount] = MatchScore;
            SecondCount++;
        }
        FirstCount++;
    }
    return Score;
}


NJOI::NJOI (std::istream* pios, int & threads):num_taxa_(0), num_char_(0), nthreads_(threads) {
    Sequence seq;
    std::string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    int seqcount = 0;
    // some error checking. should be in general seq reader class
    bool first = true;
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        sequences_[seq.get_id()] = seq.get_sequence();
        if (!first) {
            if (static_cast<int>(seq.get_length()) != num_char_) {
                std::cerr << "Error: sequence " << seq.get_id() << " has "
                    << seq.get_length() << " characters, was expecting " 
                    << num_char_ << "." << std::endl << "Exiting." << std::endl;
                exit(1);
            }
        } else {
            num_char_ = seq.get_length();
            first = false;
        }
        seqcount++;
    }
    //fasta has a trailing one
    if (ft == 2) {
        sequences_[seq.get_id()] = seq.get_sequence();
        if (static_cast<int>(seq.get_length()) != num_char_) {
            std::cerr << "Error: sequence " << seq.get_id() << " has "
                << seq.get_length() << " characters, was expecting " 
                << num_char_ << "." << std::endl << "Exiting." << std::endl;
            exit(1);
        };
        seqcount++;
    }
    num_taxa_ = seqcount;
    set_name_key ();
    Matrix_ = BuildMatrix(sequences_);
    TREEMAKE(names_, name_key_, Matrix_);
}


void NJOI::set_name_key () {
    int count = 0;
    for (iter_ = sequences_.begin(); iter_ != sequences_.end(); iter_++) {
        name_key_[count] = iter_ -> first;
        names_.push_back(iter_ -> first);
        count++;
    }
}


std::string NJOI::get_newick () {
    return newick_string_;
}
