#include <string>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
#include <map>
#include <iterator>
#include <iomanip>
#include <iostream>

using namespace std;

#include "ls_Seq.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

//Much less dumb than the way I did the nucleotides
void Stats::AA_STAT_Getter(string& seq){

	AA_Total.clear();
	string AA_string = "ACDEFGHIKLMNPQRSTVWXY*";
	for (int i = 0; i < AA_string.length(); i++){
			AA_Total[AA_string[i]] = 0.0;
	}
	for (int i = 0; i < seq.length(); i++){
		seq[i] = toupper(seq[i]);
		AA_Total[seq[i]]++;
	}
}
void Stats::GC_Getter(string& seq){

	G = 0.0; 
	C = 0.0;
	A = 0.0;
	T = 0.0;
	Missing = 0.0;
	for (int i = 0; i < seq.length(); i++){
		
			if (seq[i] == 'A' or seq[i] == 'a'){
				
					A++;
			}else if (seq[i] == 'C' or seq[i] == 'c'){
				
					C++;
			}else if (seq[i] == 'G' or seq[i] == 'g'){
				
					G++;
			}else if (seq[i] == 'T' or seq[i] == 't'){
				
					T++;
			}else{
				
				Missing++;
			}
	}
		
}

Stats::Stats (istream* pios, bool& all, bool& prot) {

	//Concatenated will be used for all stats
	const char separator    = ' ';
	const int nameWidth     = 10;
	const int numWidth      = 10;
	if (prot == false){
		Sequence seq;
		string retstring;
		int ft = test_seq_filetype_stream(*pios, retstring);
		int seqcount = 0;
		bool first = true;
		string type = "";
		while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
		
			seqcount++;
			Concatenated += seq.get_sequence();
			if (all == true){
			
				cout << "Stats For " << seq.get_id() << endl;
				temp_seq = seq.get_sequence();
				cout << "Length of Sequence: " << temp_seq.length() << endl;
				GC_Getter(temp_seq);
				GC_Content = G + C;
				Percent_A = A / temp_seq.length();
				Percent_C = C / temp_seq.length();
				Percent_G = G / temp_seq.length();
				Percent_T = T / temp_seq.length();
				Percent_GC = GC_Content / temp_seq.length();
				Percent_Missing = Missing / temp_seq.length();
				cout << "----------------------Nucleotide Composition Table-------------------------" << endl;
				cout << left << setw(nameWidth) << setfill(separator) << "\t    A";
				cout << left << setw(nameWidth) << setfill(separator) << "    C";
				cout << left << setw(numWidth) << setfill(separator) << "    G";
				cout << left << setw(numWidth) << setfill(separator) << "    T";
				cout << left << setw(numWidth) << setfill(separator) << "    G+C";
				cout << left << setw(numWidth) << setfill(separator) << "   Ambiguous";
				cout << endl;
				cout << left << setw(nameWidth) << setfill(separator) << "Total: ";
				cout << left << setw(nameWidth) << setfill(separator) << A;
				cout << left << setw(nameWidth) << setfill(separator) << C;
				cout << left << setw(numWidth) << setfill(separator) << G;
				cout << left << setw(numWidth) << setfill(separator) << T;
				cout << left << setw(numWidth) << setfill(separator) << GC_Content;
				cout << left << setw(numWidth) << setfill(separator) << Missing;
				cout << endl;
				cout << left << setw(nameWidth) << setfill(separator) << "Percent: ";
				cout << left << setw(nameWidth) << setfill(separator) << Percent_A;
				cout << left << setw(nameWidth) << setfill(separator) << Percent_C;
				cout << left << setw(numWidth) << setfill(separator) << Percent_G;
				cout << left << setw(numWidth) << setfill(separator) << Percent_T;
				cout << left << setw(numWidth) << setfill(separator) << Percent_GC;
				cout << left << setw(numWidth) << setfill(separator) << Percent_Missing;
				cout << endl;
				cout << "---------------------------------------------------------------------------" << endl;
	
			}
			if (ft == 1){
					type = "Phylip";
			}
			if (ft == 0){
					type = "Nexus";
			}
		
		}
		if (ft == 2) {
			seqcount++;
			Concatenated += seq.get_sequence();
			type = "Fasta";
			if (all == true){
			
				cout << "Stats For " << seq.get_id() << endl;
				temp_seq = seq.get_sequence();
				cout << "Length of Sequence: " << temp_seq.length() << endl;
				GC_Getter(temp_seq);
				GC_Content = G + C;
				Percent_A = A / temp_seq.length();
				Percent_C = C / temp_seq.length();
				Percent_G = G / temp_seq.length();
				Percent_T = T / temp_seq.length();
				Percent_GC = GC_Content / temp_seq.length();
				Percent_Missing = Missing / temp_seq.length();
				cout << "----------------------Nucleotide Composition Table-------------------------" << endl;
				cout << left << setw(nameWidth) << setfill(separator) << "\t    A";
				cout << left << setw(nameWidth) << setfill(separator) << "    C";
				cout << left << setw(numWidth) << setfill(separator) << "    G";
				cout << left << setw(numWidth) << setfill(separator) << "    T";
				cout << left << setw(numWidth) << setfill(separator) << "    G+C";
				cout << left << setw(numWidth) << setfill(separator) << "   Ambiguous";
				cout << endl;
				cout << left << setw(nameWidth) << setfill(separator) << "Total: ";
				cout << left << setw(nameWidth) << setfill(separator) << A;
				cout << left << setw(nameWidth) << setfill(separator) << C;
				cout << left << setw(numWidth) << setfill(separator) << G;
				cout << left << setw(numWidth) << setfill(separator) << T;
				cout << left << setw(numWidth) << setfill(separator) << GC_Content;
				cout << left << setw(numWidth) << setfill(separator) << Missing;
				cout << endl;
				cout << left << setw(nameWidth) << setfill(separator) << "Percent: ";
				cout << left << setw(nameWidth) << setfill(separator) << Percent_A;
				cout << left << setw(nameWidth) << setfill(separator) << Percent_C;
				cout << left << setw(numWidth) << setfill(separator) << Percent_G;
				cout << left << setw(numWidth) << setfill(separator) << Percent_T;
				cout << left << setw(numWidth) << setfill(separator) << Percent_GC;
				cout << left << setw(numWidth) << setfill(separator) << Percent_Missing;
				cout << endl;
				cout << "---------------------------------------------------------------------------" << endl;
	
			}
		}
		GC_Getter(Concatenated);
		GC_Content = G + C;
		Percent_A = A / Concatenated.length();
		Percent_C = C / Concatenated.length();
		Percent_G = G / Concatenated.length();
		Percent_T = T / Concatenated.length();
		Percent_GC = GC_Content / Concatenated.length();
		Percent_Missing = Missing / Concatenated.length();
	
		cout << "General Stats For All Sequences" << endl;
		cout << "File Type: " << type << endl;
		cout << "Number of Sequences: " << seqcount << endl;
		cout << "Total Length of All Combined: " << Concatenated.length() << endl;
		cout << "----------------------Nucleotide Composition Table-------------------------" << endl;
		cout << left << setw(nameWidth) << setfill(separator) << "\t    A";
		cout << left << setw(nameWidth) << setfill(separator) << "    C";
		cout << left << setw(numWidth) << setfill(separator) << "    G";
		cout << left << setw(numWidth) << setfill(separator) << "    T";
		cout << left << setw(numWidth) << setfill(separator) << "    G+C";
		cout << left << setw(numWidth) << setfill(separator) << "   Ambiguous";
		cout << left << setw(numWidth) << setfill(separator) << "   All";
		cout << endl;
		cout << left << setw(nameWidth) << setfill(separator) << "Total: ";
		cout << left << setw(nameWidth) << setfill(separator) << A;
		cout << left << setw(nameWidth) << setfill(separator) << C;
		cout << left << setw(numWidth) << setfill(separator) << G;
		cout << left << setw(numWidth) << setfill(separator) << T;
		cout << left << setw(numWidth) << setfill(separator) << GC_Content;
		cout << left << setw(numWidth) << setfill(separator) << Missing;
		cout << left << setw(numWidth) << setfill(separator) << Concatenated.length();
		cout << endl;
		cout << left << setw(nameWidth) << setfill(separator) << "Percent: ";
		cout << left << setw(nameWidth) << setfill(separator) << Percent_A;
		cout << left << setw(nameWidth) << setfill(separator) << Percent_C;
		cout << left << setw(numWidth) << setfill(separator) << Percent_G;
		cout << left << setw(numWidth) << setfill(separator) << Percent_T;
		cout << left << setw(numWidth) << setfill(separator) << Percent_GC;
		cout << left << setw(numWidth) << setfill(separator) << Percent_Missing;
		cout << left << setw(numWidth) << setfill(separator) << "100.0";
		cout << endl;
		cout << "---------------------------------------------------------------------------" << endl;
	}else{
		//This is what to do if it is an amino acid sequence
		Sequence seq;
		string retstring;
		int ft = test_seq_filetype_stream(*pios, retstring);
		int seqcount = 0;
		bool first = true;
		string type = "";
		while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
			
			seqcount++;
			Concatenated += seq.get_sequence();
			if (all == true){
				temp_seq = seq.get_sequence();
				AA_STAT_Getter(temp_seq);
				string AA_string = "ACDEFGHIKLMNPQRSTVWXY*";
				cout << "Stats For " << seq.get_id() << endl;
				cout << "Length of Sequence: " << temp_seq.length() << endl;
				cout << "------------AA TABLE------------" << endl;
				cout << "AA\tTotal\tPercent" << endl;
				for (int i = 0; i < AA_string.length(); i++){
					cout << left << setw(nameWidth) << setfill(separator) << AA_string[i] << AA_Total[AA_string[i]] << "\t" << ((AA_Total[AA_string[i]] / temp_seq.length())*100.0) << endl;
				
				}
				cout << "------------AA TABLE------------" << endl;
			}
			if (ft == 1){
					type = "Phylip";
			}
			if (ft == 0){
					type = "Nexus";
			}
		
		}
		if (ft == 2) {
			seqcount++;
			Concatenated += seq.get_sequence();
			type = "Fasta";
			if (all == true){
				temp_seq = seq.get_sequence();
				AA_STAT_Getter(temp_seq);
				string AA_string = "ACDEFGHIKLMNPQRSTVWXY*";
				cout << "Stats For " << seq.get_id() << endl;
				cout << "Length of Sequence: " << temp_seq.length() << endl;
				cout << "------------AA TABLE------------" << endl;
				cout << "AA\tTotal\tPercent" << endl;
				for (int i = 0; i < AA_string.length(); i++){
					cout << left << setw(nameWidth) << setfill(separator) << AA_string[i] << AA_Total[AA_string[i]] << "\t" << ((AA_Total[AA_string[i]] / temp_seq.length())*100.0) << endl;
				
				}
				cout << "------------AA TABLE------------" << endl;
			}
		}
		cout << "General Stats For All Sequences" << endl;
		cout << "File Type: " << type << endl;
		cout << "Number of Sequences: " << seqcount << endl;
		cout << "Total Length of All Combined: " << Concatenated.length() << endl;
		AA_STAT_Getter(Concatenated);
		string AA_string = "ACDEFGHIKLMNPQRSTVWXY*";
		cout << "------------AA TABLE------------" << endl;
		cout << "AA\tTotal\tPercent" << endl;
		for (int i = 0; i < AA_string.length(); i++){
			cout << left << setw(nameWidth) << setfill(separator) << AA_string[i] << AA_Total[AA_string[i]] << "\t" << ((AA_Total[AA_string[i]] / Concatenated.length())*100.0) << endl;
		}
	}
}	
