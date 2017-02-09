/*
 * distmatrix.cpp
 *
 *  Created on: Jun 3, 2015
 *      Author: joe
 */


// NOTE: this file is not presently used

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <map>
#include <iterator>

using namespace std;

//Calculate the difference between two strings
float CalcSeqDiffs(string& sequence1, string& sequence2) {
    float score = 0;
    for (unsigned int i = 0; i < sequence1.size(); i++) {
        if (sequence1[i] != sequence2[i]) {
            score++;
        }
    }
    return score;
}

map<string,int> BuildMatrix (map <string, string>& sequences) {
    vector<string> SequenceName;
    map<string, int> Matrix;
    map <string, string>::iterator iter;
    map <string, string>::iterator iter2;
    string fasta, SeqName, MatchName;
    int count = 0;
    int FirstCount = 0;
    float MatchScore;
    vector< vector<int> > Score;
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {
        vector<int> row;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++) {
            row.push_back(0);
        }
        Score.push_back(row);
    }
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {
        fasta = iter -> second;
        SeqName = iter -> first;
        SequenceName.push_back(SeqName);
        float count = 0;
        int SecondCount = 0;
        for (iter2 = sequences.begin(); iter2 != sequences.end(); iter2++) {
            MatchScore = CalcSeqDiffs(fasta, iter2 -> second);
            MatchName = SeqName + "," + iter2 -> first;
            Matrix[MatchName] = MatchScore;
            Score[FirstCount][SecondCount] = MatchScore;
            SecondCount++;
        }
        FirstCount++;

    }
    cout << "\t";
    for (int i = 0; i < SequenceName.size(); i++) {

        cout << SequenceName[i] << "\t";
    }
    cout << endl;
    for (int i = 0; i < Score.size(); i++) {
        cout << SequenceName[i] << "\t";
        for (int j = 0; j < Score[i].size(); j++) {
            cout << Score[i][j] << "\t";
        }
        cout << endl;
    }
    return Matrix;
}

//Changes the input fasta to one line and returns it
//as a map
map <string, string> FastaToOneLine (string& fasta) {

    map <string, string> sequences;
    string line, dna, name_hold;
    ifstream readline;
    float count = 0;
    readline.open(fasta.c_str());
    if (readline.is_open()) {
        while (getline (readline, line)) {
            if (line[0] == '>') {
                if (count != 0) {
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                    sequences[name_hold] = dna;
                    dna = "";
                } else {
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                }
                name_hold = line;
            } else {
                count = 1;
                dna += line;
            }
        }
    }
    readline.close();
    sequences[name_hold] = dna;
    return sequences;
}

int main() {
    map <string, string> sequences;
    map <string, int>::iterator iter;
    string fasta, Tree;
    map<string, int> Matrix;
    cout << "Calculating Distances" << endl;
    //cin >> fasta;
    fasta = ("TestFiles/Real_Test.fa");
    sequences = FastaToOneLine(fasta);
    Matrix = BuildMatrix(sequences);
}
