/*
 * RemoveTaxa.cpp
 *
 *  Created on: Jan 7, 2015
 *      Author: joe
 */

/*
 * RemoveTaxa.cpp
 *
 *  Created on: Jan 5, 2015
 *      Author: joe
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>

using namespace std;

//Makes a vector with fasta sequences on one line
vector <string> ReadFasta(string fasta) {
    string line, dna;
    ifstream readline;
    int count = 0;
    vector<string> parsed_fasta;
    readline.open(fasta.c_str());
    if (readline.is_open()) {
        while (getline (readline, line)) {
            if (line[0] == '>') {

                if (count != 0) {
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                    parsed_fasta.push_back(dna);
                    parsed_fasta.push_back(line);
                    dna = "";
                } else {
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                    parsed_fasta.push_back(line);
                }
            } else {
                count++;
                dna += line;
            }
        }
    }
    parsed_fasta.push_back(dna);
    return parsed_fasta;
}

//Reads sequences that need to be removed
vector <string> ReadSeqs(string seqs) {
    string line, dna;
    ifstream readline;
    vector<string> parsed_seqs;
    readline.open(seqs.c_str());
    if (readline.is_open()) {
        while (getline (readline, line)) {
            line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
            parsed_seqs.push_back(line);
        }
    }
    return parsed_seqs;
}

int main() {
    string fasta;
    string ListofSeqs;
    vector<string> parsed_fasta;
    vector<string> parsed_seqs;
    cout << "Type in Name of Seq" << endl;
    cin >> fasta;
    cout << "Input a space delimited file of Sequences to be removed" << endl;
    cin >> ListofSeqs;
    parsed_fasta = ReadFasta(fasta);
    parsed_seqs = ReadSeqs(ListofSeqs);
    //checks for seqs that need to be removed and removes them
    for (vector<string>::size_type i = 0; i != parsed_fasta.size(); i = i + 2) {
        for (vector<string>::size_type j = 0; j != parsed_seqs.size(); j++) {
            if (parsed_seqs[j] == parsed_fasta[i]) {
                //cout << ">" << parsed_fasta[i] << "\n" << parsed_fasta[i+1] << endl;
                parsed_fasta.erase(remove(parsed_fasta.begin(), parsed_fasta.end(), parsed_fasta[i]), parsed_fasta.end());
                parsed_fasta.erase(remove(parsed_fasta.begin(), parsed_fasta.end(), parsed_fasta[i]), parsed_fasta.end());
            }
        }
    }
    //prints out remaining seqs
    for (vector<string>::size_type i = 0; i != parsed_fasta.size(); i = i + 2) {
        cout << ">" << parsed_fasta[i] << "\n" << parsed_fasta[i+1] << endl;
    }

    return EXIT_SUCCESS;
}
