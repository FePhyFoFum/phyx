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
#include <map>
#include <iterator>

using namespace std;

#include "njoi.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

/*Calculates the Q matrix
 * Original Matrix: All the Original Distances
 * ConvertedMatrix: Conversion to Q Matrix
 * LengthMatrix: Adjusted values for branch length calculations
 */
void CalcQ(int NumbOfSequences, vector< vector<double> >& OriginalMatrix,vector< vector<double> >& ConvertedMatrix, vector< vector<double> >& LengthMatrix){


	ConvertedMatrix = OriginalMatrix;
	double row_sum = 0.0;
	double col_sum = 0.0;
	vector<double> Sums;
	for (int i = 0; i < NumbOfSequences; i++){
		row_sum = 0.0;
		for (int j = 0; j < NumbOfSequences; j++){
				row_sum += ConvertedMatrix[i][j];
				ConvertedMatrix[i][j] *= (NumbOfSequences - 2);

		}
		Sums.push_back(row_sum);
		if (i == 0){
			col_sum = row_sum;
		}
	}
	Sums.push_back(col_sum);
	for (int i = 0; i < NumbOfSequences; i++){
		for (int j = 0; j < NumbOfSequences; j++){

			if (i != j){
				LengthMatrix[i][j] = abs(Sums[i] - Sums[j]);
				ConvertedMatrix[i][j] -= Sums[i];
				ConvertedMatrix[i][j] -= Sums[j];
			}


		}
	}
	Sums.clear();
}
/*Calculate the Branch Lengths
 * NewMatrix: Has original distances
 * LengthMatrix: Has adjusted lengths
 */
void FetchLengths(int NumbOfSequences, vector< vector<double> >& NewMatrix, vector< vector<double> >& LengthMatrix,
		int& mini1, int& mini2, double& brlength1, double& brlength2){

	//Cant figure out how to bypass changing to a string then a double, but their has to
	//be a command to change integer to double
	string hold = to_string(NumbOfSequences);
    brlength1 = (NewMatrix[mini1][mini2] / 2.0) + (1 / (2*(stod(hold)-2))) * LengthMatrix[mini1][mini2];
    brlength2 = NewMatrix[mini1][mini2] - brlength1;

}
/*Tree_Update
 *Updates the Tree info and bypasses using a tree structure by storing parts in an array
 *NewMatrix: Has the adjusted values from the QMatrix calculation
 */

void Tree_Update(int NumbOfSequences, string& newname, vector<string>& names, map<int, string>& NumbKeys,
    int& node_list, vector< vector<double> >& NewMatrix, int& mini1, int& mini2, double& brlength1, double& brlength2){

    //update the tree values, Tree Size is the node it is at
    vector<double> row_hits, col_hits, new_ColRow;
    double ColRow = 0.0;
    double small_length = NewMatrix[mini1][mini2]; // neighbor based correction
    newname = "(" + names[mini1] + ":" + to_string(brlength2) +  "," + names[mini2] + ":" + to_string(brlength1) +  ")";
    names.erase(names.begin()+mini1);
    names.erase(names.begin()+(mini2-1));
    names.insert(names.begin(), newname);
    //Make Smaller Matrix
    vector< vector<double> > temp_matrix(node_list, vector<double>(node_list, 0.0));
    //Reformat Matrix Takes in the hits
    for (int i = 0; i < NewMatrix.size(); i++){
            for (int j = 0; j < NewMatrix.size(); j++){
                if (i == mini1){
                        row_hits.push_back(NewMatrix[mini1][j]);
                }else if(i == mini2){
                    col_hits.push_back(NewMatrix[mini2][j]);
                }
            }
    }
    int count = 0;
    //Make a new First Row and Column
    for (int i = 0; i < col_hits.size(); i++){
        ColRow = (col_hits[i] + row_hits[i] - small_length) / 2;
        if (i != mini1){
            if (i != mini2){
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
    for (int i = 0; i < NewMatrix.size(); i++){
        jcount = 1;
        if (i != mini1){
            if(i != mini2){
                for (int j = 0; j < NewMatrix.size(); j++){
                    if (j != mini1){
                        if(j != mini2){
                            temp_matrix[icount][jcount] = NewMatrix[i][j];
                            jcount++;
                        }
                    }
                }
                icount++;
            }
        }
    }
    NewMatrix = temp_matrix;

}

void Choose_Smallest(int& node_list, vector< vector<double> >& Matrix, int& mini1, int& mini2){

    //super large value
    double MIN = 99999999999.99;
    for (int i = 0; i < node_list; i++){
        for(int j = i + 1; j < node_list; j++){
            if (Matrix[i][j] < MIN){
                MIN = Matrix[i][j];
                mini1 = i;
                mini2 = j;
            }
        }
    }
    node_list--;
}

//NumbKeys Contains the names and their matching number
//Matrix Contains The original matrix
//Main Tree Making Matrix
void NJOI::TREEMAKE(vector<string>& names, map <int, string>& NumbKeys,
    vector< vector<double> >& Matrix){

    int mini1 = 0, mini2 = 0;
    int NumbOfSequences = NumbKeys.size();
    double brlength1 = 0.0;
    double brlength2 = 0.0;
	vector< vector<double> > NewMatrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
	vector< vector<double> > LengthMatrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
    NewMatrix = Matrix;
    map<int, string>::iterator iter;
    string newname;
    while (NumbOfSequences > 2){


    	vector< vector<double> > QMatrix(NumbOfSequences, vector<double>(NumbOfSequences, 0.0));
        CalcQ(NumbOfSequences, Matrix, QMatrix, LengthMatrix);
        Choose_Smallest(NumbOfSequences, QMatrix, mini1, mini2);
        FetchLengths((NumbOfSequences + 1), Matrix, LengthMatrix, mini1, mini2, brlength1, brlength2);
        Tree_Update((NumbOfSequences + 1), newname, names, NumbKeys, NumbOfSequences, Matrix, mini1, mini2, brlength1, brlength2);

    }
    double adjlength = (Matrix[mini1][mini2] / 2); // The final branch length
    //cout << Matrix[mini1][mini2] << endl;
    //FetchLengths((NumbOfSequences + 1), Matrix, LengthMatrix, mini1, mini2, brlength1, brlength2);
    newname = "(" + names[mini1] + ":" + to_string(adjlength) +  "," + names[mini2] + ":" + to_string(adjlength) +  ")";
    newickstring = newname + ";";
}



double CalcSeqsDiffs(string& sequence1, string& sequence2){

    double score = 0;
    for (int i = 0; i < sequence1.size(); i++){

        if (sequence1[i] != sequence2[i]){

            score++;
        }
    }
    return score;
}

vector< vector<double> > NJOI::BuildMatrix (map <string, string>& sequences){

    vector<string> SequenceName;
    map <string, string>::iterator iter, iter2;
    string fasta, SeqName, MatchName;
    int count = 0, FirstCount = 0;
    double MatchScore;

    // an easier way to initialize a vector of vectors:
    int ntax = sequences.size();
    vector< vector<double> > Score(ntax, vector<double>(ntax, 0.0));

    //compare all sequences to other sequences
    for (iter = sequences.begin(); iter != sequences.end(); iter++){
        fasta = iter -> second;
        SeqName = iter -> first;
        SequenceName.push_back(SeqName);
        int SecondCount = 0;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++){

            MatchScore = CalcSeqsDiffs(fasta, iter2 -> second);
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

NJOI::NJOI (string & seqf):ntax(0), nchar(0) {
    ifstream* fstr;
    istream* pios;
    fstr = new ifstream(seqf);
    pios = fstr;
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios,retstring);
    while(read_next_seq_from_stream(*pios,ft,retstring,seq)){
        sequences[seq.get_id()] = seq.get_sequence();
    }
    //fasta has a trailing one
    if (ft == 2){
    	sequences[seq.get_id()] = seq.get_sequence();
    }
    ntax = sequences.size();
    nchar = (int)sequences.begin()->second.size(); // ugly. should set when reading seqs
    set_name_key ();
    Matrix = BuildMatrix(sequences);
    TREEMAKE(names, NameKey, Matrix);
}

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
