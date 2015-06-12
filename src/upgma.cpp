/*
 * upgma.cpp
 *
 *  Created on: Jun 10, 2015
 *      Author: joe
 */

//#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <map>
#include <iterator>

using namespace std;

#include "upgma.h"
#include "utils.h"

void Update_Tree(string& newname, vector<string>& names, map<int, string>& NumbKeys, 
    int& node_list, vector< vector<double> >& NewMatrix, int& mini1, int& mini2){

    //update the tree values, Tree Size is the node it is at
    vector<double> row_hits, col_hits, new_ColRow;
    double br_length = NewMatrix[mini1][mini2] / 2.0;
        
    string length1 = std::to_string(br_length);
    string length2 = length1;
    double ColRow = 0.0;

    // extremely hacky way to get correct ELs
    std::size_t found = names[mini1].find("#");
    if (found != std::string::npos) {
        string terp = names[mini1];
        std::size_t pos = terp.find("#"); 
        double oldheight = std::stod(terp.substr(pos+1));
        names[mini1] = terp.substr(0, pos);
        length1 = std::to_string(br_length - oldheight);
    }
    // have to do it for the other side too
    found = names[mini2].find("#");
    if (found != std::string::npos) {
        string terp = names[mini2];
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
    vector< vector<double> > temp_matrix;
        
    for (int i = 0; i < node_list; i++){
        vector<double> row;
        for (int j = 0; j < node_list; j++){
            row.push_back(0);
                }
        temp_matrix.push_back(row);
    }
    //Reformat Matrix
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
        ColRow = (col_hits[i] + row_hits[i]) / 2;
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
    //Print the new matrix makes it very verbose but
    //fun to watch
    /*
    cout << name1 << "  <===============NewMatrix==============>  " << name2 << endl;
    for (int i = 0; i < temp_matrix.size(); i++){
        for (int j = 0; j < temp_matrix.size(); j++){
            cout <<temp_matrix[i][j] << "\t";
        }
        cout << endl;
    }*/
    NewMatrix = temp_matrix;

}

void Choose_Small(int& node_list, vector< vector<double> >& Matrix, int& mini1, int& mini2){

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
void UPGMA::TREEMAKE(vector<string>& names, map <int, string>& NumbKeys,
    vector< vector<double> >& Matrix){

    int mini1 = 0, mini2 = 0;
    int NumbOfSequences = NumbKeys.size();
    vector< vector<double> > NewMatrix;
    map<int, string>::iterator iter;
    string newname;
    while (NumbOfSequences > 1){

        Choose_Small(NumbOfSequences, Matrix, mini1, mini2);
        Update_Tree(newname, names, NumbKeys, NumbOfSequences, Matrix, mini1, mini2);

    }
    newickstring = newname + ";";
}



double CalcSeqDiffs(string& sequence1, string& sequence2){

    double score = 0;
    for (int i = 0; i < sequence1.size(); i++){

        if (sequence1[i] != sequence2[i]){

            score++;
        }
    }
    return score;
}

vector< vector<double> > UPGMA::BuildMatrix (map <string, string>& sequences){

    vector<string> SequenceName;
    map <string, string>::iterator iter, iter2;
    string fasta, SeqName, MatchName;
    int count = 0, FirstCount = 0;
    double MatchScore;
        
    //vector< vector<double> > Score;
    //initialize the 2d vector and fill it with zeroes
        /*
    for (iter = sequences.begin(); iter != sequences.end(); iter++){
        vector<double> row;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++){
            row.push_back(0);
                }
        Score.push_back(row);
    }
        */
       
    // an easier way to initialize a vector of vectors:
    int ntax = sequences.size(); // this should be a property of the class, calculated once
    vector< vector<double> > Score(ntax, vector<double>(ntax, 0.0));
    
    //compare all sequences to other sequences
    for (iter = sequences.begin(); iter != sequences.end(); iter++){
        fasta = iter -> second;
        SeqName = iter -> first;
        SequenceName.push_back(SeqName);
        int SecondCount = 0;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++){

            MatchScore = CalcSeqDiffs(fasta, iter2 -> second);
            MatchName = SeqName + "," + iter2 -> first;
            Score[FirstCount][SecondCount] = MatchScore;
            SecondCount++;

        }
        FirstCount++;

    }
    //prints the distance matrix maybe too verbose
    /*
    cout << "\t";
    for (int i = 0; i < SequenceName.size(); i++){

        cout << SequenceName[i] << "\t";
    }
    cout << endl;
    for(int i = 0; i < Score.size(); i++){
        cout << SequenceName[i] << "\t";
        for (int j = 0; j < Score[i].size(); j++){
            cout << Score[i][j] << "\t";

        }
        cout << endl;
    }*/
    return Score;
}

// should use existing seq readers
map<string, string> UPGMA::FastaToOneLine (string& fasta){

    string line, dna, name_hold;
    ifstream readline;
    bool round_one = true;
    readline.open(fasta.c_str());
    if (readline.is_open()){
        while (getline (readline, line)){
            if (line[0] == '>'){

                if (round_one == false){
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                    sequences[name_hold] = dna;
                    dna = "";

                }else{
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                }
                name_hold = line;
            }else{
                round_one = false;
                dna += line;
            }
        }
    }
    sequences[name_hold] = dna;
    return sequences;

}

UPGMA::UPGMA() {
    // TODO Auto-generated constructor stub

}

// *** some alternate functions below *** //

// alternate constructor
UPGMA::UPGMA (string & seqf):ntax(0), nchar(0) {
    seqfile = seqf;
    sequences = FastaToOneLine (seqfile); // use existing seq readers instead
    ntax = sequences.size();
    nchar = (int)sequences.begin()->second.size(); // ugly. should set when reading seqs
    set_name_key ();
    Matrix = BuildMatrix(sequences);
    TREEMAKE(names, NameKey, Matrix);
}

// just moving this over from main
void UPGMA::set_name_key () {
    int count = 0;
    for(iter = sequences.begin(); iter != sequences.end(); iter++) {
        NameKey[count] = iter -> first;
        names.push_back(iter -> first);
        count++;
    }
}

string UPGMA::get_newick () {
    return newickstring;
}