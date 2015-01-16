/*
 * RemoveAmbiguous.cpp
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
//Get number of lines
int numb_lines(string dna){

	int count = 0;
	string name, seq, dna_file;
	ifstream readline;
	vector<int> locations;

	readline.open(dna.c_str());
	if (readline.is_open()){

		while(readline >> name >> seq){

			count++;
		}
	}
	return (count - 1);

}

//Parse file and where N appears in each sequence

vector<int> count_ambiguous(string dna){

	int count = 0;
	string name, seq, dna_file, line;
	ifstream readline;
	vector<int> locations;

	readline.open(dna.c_str());
	if (readline.is_open()){

		while(readline >> name >> seq){

			for (int i = 0; i < seq.size(); i++){
			    if (seq[i] == 'N'){

			    	locations.push_back(i);
			    }
			}
		}
	}
	return locations;


}
//take the vector with N's, sort it and find how many times each number appears
vector<int> sorter(vector<int>& dna){

	vector<int> sorted_vector;
	vector<int> return_vec;
	sort( dna.begin(), dna.end() );

	for (vector<int>::const_iterator it=dna.begin(); it!=dna.end(); ++it) {
	     // cout << *it << " ";
	      sorted_vector.push_back(*it);
	}
	for(size_t i = 0; i < sorted_vector.size(); i++)
	{
	    size_t count = 1;

	    size_t limit = sorted_vector.size() - 1;
	    while(i < limit  && sorted_vector[i] == sorted_vector[i+1])
	    {
	    	count++;
	    	i++;
	    }
	    return_vec.push_back(sorted_vector[i]);
	    return_vec.push_back(count);
	}
	return return_vec;

}

//Find positions greater than cutoff
vector<int> cutoff_identifier(vector<int>& sorted_vector, float cutoff, float number_of_seqs){

	vector<int> problem_positions;
	float percent_occurance = 0;
	for(vector<string>::size_type i = 0; i != sorted_vector.size(); i = i + 2) {

		//cout << "Postion: " << sorted_vector[i] << "\tOcurrances: " << sorted_vector[i + 1] << '\n';
		percent_occurance = sorted_vector[i + 1] / number_of_seqs;
		//cout << "\t" << sorted_vector[i+ 1] << "\t" << number_of_seqs << "\t" << percent_occurance << endl;
		if (cutoff < percent_occurance){
			//cout << sorted_vector[i] << "\t" << sorted_vector[i + 1] << endl;
			problem_positions.push_back(sorted_vector[i]);

		}

	}
	return problem_positions;
}

string trimmer(vector<int>& problem_positions, string seq){

	string cat;
	for (int i = 0; i < seq.size(); i++){
		if(std::find(problem_positions.begin(), problem_positions.end(), i) != problem_positions.end()) {
	    	/* problem positions contains i */
			//cout << i << ' ';
			//cout << seq[i];
		}else{
	    	cat = cat + seq[i];
		}
	}
	return cat;
}

int main()
{
	string dna_file, name, seq, trimmed_seq; //Input file
	vector<int> positions; //contains position of each N
	vector<int> sorted_vector, problem_positions;
	float cutoff = 0;
	float number_of_seqs = 0; //contains number of seqs in file
	ifstream readline;
	ofstream output;
	output.open("trimmed.phy");

	//Read in the Phylip File

	cout << "Please Enter File" << endl;
	cin >> dna_file;

	//Read in the percent Ambiguous to remove

	cout << "Percent cutoff for ambiguous sequences" << endl;
	cin >> cutoff;
	cout << "If more than " << cutoff << "% of your bases at a position"
		 << " are ambiguous that position will be removed" << endl;
	cutoff = cutoff / 100;
	//Get positions of N's
	positions = count_ambiguous(dna_file);
	//Get number of seqs
	number_of_seqs = numb_lines(dna_file);
	//take the vector with N's and find how many times each number appears
	sorted_vector = sorter(positions);
	//Find the positions greater than the cutoff, stor in vector
	problem_positions = cutoff_identifier(sorted_vector, cutoff, number_of_seqs);
	//Parse file again with counts of N's and remove the N's
	readline.open(dna_file.c_str());
	if (readline.is_open()){

		while(readline >> name >> seq){

			trimmed_seq = trimmer(problem_positions, seq);

			cout << name << "\t" << trimmed_seq << endl;
			output << name << "\t" << trimmed_seq << endl;
		}
	}
	output.close();

}
