/*
 * concat.cpp
 *
 *  Created on: Sep 22, 2014
 *      Author: joe
 */



#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
//#include <string.h>
#include <algorithm>
//#include <set>
using namespace std;

//Gets the length of the string
int getlength(string seq){

	int seq_length;
	seq_length = seq.size();
	return seq_length;

}
//function splits files on tab and returns vector
//each element of the vector contains what has been
//split on the tab, so element 0 is name and
//element 1 is the sequence itself
vector<string> tab_split(string dna){

	string second_dna_file;
	int count = 0;
	int Length_of_token;
	string line;
	ifstream input_dna;
	vector<string> names; //vector to contain names
	vector<string> sequences; //vector to contain sequences
	string token;
	input_dna.open(dna.c_str());
		if (input_dna.is_open()){

			//split into separate strings
			while (getline (input_dna, line)){
				//cout << line << endl;
				istringstream str(line);
				if (count != 0){
					while(getline(str, token, ' ')){
						//create a vector with name and sequence as elements
						 //   token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
							Length_of_token = getlength(token);
							if (Length_of_token != 0){
								cout << "After Erase" << token << '\n';
								names.push_back(token);
						//	cout << token << endl;
							}
					}
				}
				count++;
			}
		}
		return names;
}

//not a good way but the best I can come up with now
//to return all the unique values of a vector
vector<string> matcher(vector<string>& unique_ones, vector<string>& match)
{
	//sorts uniquely and removes duplicates
	 vector<string>::iterator it;
	 vector<string> no_dup;
     sort(unique_ones.begin(), unique_ones.end());
     it = unique(unique_ones.begin(), unique_ones.end());
     unique_ones.resize(distance(unique_ones.begin(),it));
     unique_ones.push_back("@!#!@$!@$!@$!@");

     for(vector<string>::size_type i = 0; i != match.size(); i++) {
			for(vector<string>::size_type j = 0; j != unique_ones.size(); j++){
				if(unique_ones[j] == match[i]){
					unique_ones.erase(remove(unique_ones.begin(), unique_ones.end(), unique_ones[j]), unique_ones.end());
				}

			}
		}
     return unique_ones;
}

//where the strings get concatenated
//not great but doesn't suck nearly bad as attempts 1 or 2
//vector inputs are the running string called sequences and
//the new string called second_sequences
//integer inputs are the running length of a string and
//the length of the newest string
vector<string> concat(vector<string>& old_vec, vector<string>& new_vec, int current_size, int running_size){

	vector<string> phylip;
	vector<string> no_match_old;
	vector<string> no_match_new;
	vector<string> match;
	vector<string> unique_old;
	vector<string> unique_new;
	int old_match_count = 0;
	int new_match_count = 0;
	string cat;
	for(vector<string>::size_type i = 0; i != old_vec.size(); i = i + 2) {

		cout << "old vector: " << old_vec[i] << "\t" << old_vec[i + 1] << '\n';
	}
	for(vector<string>::size_type i = 0; i != new_vec.size(); i = i + 2) {

		cout << "new vector: " << new_vec[i] << "\t" << new_vec[i + 1] << '\n';
	}

	//find matches of old vector to new
	for(vector<string>::size_type i = 0; i != old_vec.size(); i = i + 2) {
			for(vector<string>::size_type j = 0; j != new_vec.size(); j = j + 2){
				if (old_vec[i] == new_vec[j]){
					cat = old_vec[i + 1] + new_vec[j+1];
					phylip.push_back(old_vec[i]);
					phylip.push_back(cat);
					match.push_back(old_vec[i]);

				}else{
					//filled to the top with non_matches
					no_match_old.push_back(old_vec[i]);
					no_match_new.push_back(new_vec[j]);
				}
			}
	}
	for(vector<string>::size_type i = 0; i != old_vec.size(); i = i + 2) {
		for(vector<string>::size_type j = 0; j != match.size(); j++){
			if (match[j] == old_vec[i]){

				old_match_count++;
			}
		}

	}
	for(vector<string>::size_type i = 0; i != new_vec.size(); i = i + 2) {
		for(vector<string>::size_type j = 0; j != match.size(); j++){
			if (match[j] == new_vec[i]){

				new_match_count++;
			}
		}

	}
	if ((old_vec.size() / 2) != old_match_count){
		unique_old = matcher(no_match_old, match);
		for(vector<string>::size_type i = 0; i != unique_old.size(); i++) {
			cout << "Unique old ones: " << unique_old[i] << '\n';
		}
		for(vector<string>::size_type i = 0; i != unique_old.size(); i++) {
			for(vector<string>::size_type j = 0; j != old_vec.size(); j = j + 2){
				if (unique_old[i] == old_vec[j]){
					string stuff(current_size, 'N');
					cat = old_vec[j + 1] + stuff;
					phylip.push_back(old_vec[j]);
					phylip.push_back(cat);

				}
			}
		}
	}
	if ((new_vec.size() / 2) != new_match_count){
		unique_new = matcher(no_match_new, match);
		for(vector<string>::size_type i = 0; i != unique_new.size(); i++) {
			for(vector<string>::size_type j = 0; j != new_vec.size(); j = j + 2){
				if(unique_new[i] == new_vec[j]){

					string stuff(running_size, 'N');
					cat = stuff + new_vec[j + 1];
					phylip.push_back(new_vec[j]);
					phylip.push_back(cat);

				}

			}
		}
	}
	return phylip;
}
int main()
{

	string dna;
	string second_dna_file;
	string line;
	string cat;
	char cont = 'Y';
	int while_count = 0;
	int running_size = 0;
	int current_size = 0;
	int total_size = 0;
	int count = 0;
	vector<string> together;
	ifstream input_dna;
	ofstream concatenated;
	ofstream final_product;
	concatenated.open("sequence3.txt");
	final_product.open("final.phy");
	vector<string> temporary; //vector to contain names
	vector<string> sequences; //vector to contain sequences
	vector<string> second_sequences;
	string token;
	//Here is where the file is read in line by line
	cout << "enter dna: ";
	cin >> dna;
	//have new sequences enter here

	while(cont == 'Y'){

		cout << "enter next: ";
		cin >> second_dna_file;
		cout << "continue: ";
		cin >> cont;
		if(while_count == 0){

			sequences = tab_split(dna);
			running_size = getlength(sequences[1]);
			concatenated << "The size is: " << running_size << endl;
			second_sequences = tab_split(second_dna_file);
			current_size = getlength(second_sequences[1]);
			concatenated << "The size of new sequence is: " << current_size << endl;
			together = concat(sequences, second_sequences, current_size, running_size);
			while_count++;
			total_size = current_size + running_size;
			concatenated << "The running size is: " << total_size << endl;

		}else{
			second_sequences = tab_split(second_dna_file);
			current_size = getlength(second_sequences[1]);
			together = concat(together, second_sequences, current_size, total_size);
			total_size = current_size + total_size;
			concatenated << "The size is: " << current_size << endl;
			concatenated << "The running size is: " << total_size << endl;
		}
		for(vector<string>::size_type i = 0; i != together.size(); i++){

			cout << together[i] << "\n";
		}

	}
	//Here is where the phylip file is printed out
	final_product << (together.size() / 2) << " " << total_size << "\n";
	for(vector<string>::size_type i = 0; i != together.size(); i++){

		if (count == 1){
			final_product << together[i] << "\n";
			count = 0;
		}else{
			final_product << together[i] << "\t";
			count = 1;
		}
	}
	//Here is where they are closed
	final_product.close();
	concatenated.close();
	return 0;

}





