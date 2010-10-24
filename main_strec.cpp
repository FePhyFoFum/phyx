/*
 * main_strec.cpp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "sequence.h"
#include "phylip_reader.h"
#include "state_reconstructor.h"
#include "rate_model.h"
#include "optimize_state_reconstructor_nlopt.h"

int main(int argc, char * argv[]){
	TreeReader tr;

	if (argc != 2){
		cout << "usage: phyx_stmap configfile" << endl;
		exit(0);
	}

	bool verbose = true;
	string datafile;
	string treefile;
	string outfile_anc = "";
	string outfile_stochtime = "";
	string outfile_stochnum ="";
	string outfile_stochnum_any ="";
	bool datawide = false;
	map<string,vector<string> > mrcas;
	vector<string> stochtime;
	vector<string> stochnumber;
	vector<string> stochnumber_any;
	vector<string> ancstates;
	string freeparams = "_one_"; //right now, just _all_ or _one_

	/*************
	 * read the configuration file
	 **************/
	ifstream ifs(argv[1]);
	string line;
	while(getline(ifs,line)){
		if(line.size()>1){
			if((&line[0]) != "#"){
				vector<string> tokens;
				string del("=");
				tokens.clear();
				Tokenize(line, tokens, del);
				for(unsigned int j=0;j<tokens.size();j++){
					TrimSpaces(tokens[j]);
				}
				if(!strcmp(tokens[0].c_str(), "treefile")){
					treefile = tokens[1];
				}else if(!strcmp(tokens[0].c_str(),  "datafile")){
					datafile = tokens[1];
				}else if(!strcmp(tokens[0].c_str(),  "freeparams")){
					freeparams = tokens[1];
				}else if(!strcmp(tokens[0].c_str(),  "outanc")){
					outfile_anc = tokens[1];
				}else if(!strcmp(tokens[0].c_str(),  "outsttime")){
					outfile_stochtime = tokens[1];
				}else if(!strcmp(tokens[0].c_str(),  "outstnum_any")){
					outfile_stochnum_any = tokens[1];
				}else if(!strcmp(tokens[0].c_str(),  "outstnum")){
					outfile_stochnum = tokens[1];
				}else if(!strcmp(tokens[0].c_str(),  "ratematrix")){
					ratematrixfile = tokens[1];
					if(ratematrixfile == "d" || ratematrixfile == "D"){
						ratematrixfile = "";
					}
				}else if(!strcmp(tokens[0].c_str(), "mrca")){
					vector<string> searchtokens;
					Tokenize(tokens[1], searchtokens, ", 	");
					for(unsigned int j=0;j<searchtokens.size();j++){
						TrimSpaces(searchtokens[j]);
					}
					vector<string> mrc;
					for(unsigned int j=1;j<searchtokens.size();j++){
						mrc.push_back(searchtokens[j]);
					}
					mrcas[searchtokens[0]] = mrc;
				}else if(!strcmp(tokens[0].c_str(), "stochtime")){
					vector<string> searchtokens;
					Tokenize(tokens[1], searchtokens, ", 	");
					for(unsigned int j=0;j<searchtokens.size();j++){
						TrimSpaces(searchtokens[j]);
						stochtime.push_back(searchtokens[j]);
					}
				}else if(!strcmp(tokens[0].c_str(), "stochnumber")){
					vector<string> searchtokens;
					Tokenize(tokens[1], searchtokens, ", 	");
					for(unsigned int j=0;j<searchtokens.size();j++){
						TrimSpaces(searchtokens[j]);
						stochnumber.push_back(searchtokens[j]);
					}
				}else if(!strcmp(tokens[0].c_str(), "stochnumber_any")){
					vector<string> searchtokens;
					Tokenize(tokens[1], searchtokens, ", 	");
					for(unsigned int j=0;j<searchtokens.size();j++){
						TrimSpaces(searchtokens[j]);
						stochnumber_any.push_back(searchtokens[j]);
					}
				}else if(!strcmp(tokens[0].c_str(), "ancstates")){
					vector<string> searchtokens;
					Tokenize(tokens[1], searchtokens, ", 	");
					for(unsigned int j=0;j<searchtokens.size();j++){
						TrimSpaces(searchtokens[j]);
						ancstates.push_back(searchtokens[j]);
					}
				}else if(!strcmp(tokens[0].c_str(),  "datawide")){
					datawide = true;
				}
			}
		}
	}
	if(verbose)
		cout << "finished reading config file" << endl;


	vector<Sequence> seqs;
	PhylipReader pr;
	bool phyl = pr.readFile(datafile,seqs);
	if(phyl == false){
		cout << "the sequence file is not phylip" << endl;
		exit(0);
	}
	if(verbose)
		cout << "sequences: " << seqs.size() << endl;

	/////////////////

	ifstream infile2(treefile.c_str());
	if (!infile2){
		cerr << "Could not open treefile." << endl;
		return 1;
	}
	vector<string> lines;
	line = "";
	while (getline(infile2, line)){
		if(line.length() > 5)
			lines.push_back(line);
	}
	infile2.close();

	///////////////

	int nstates;
	int nsites;
	if (datawide){
		nstates = seqs[0].get_sequence().length();
		nsites = 1;
	}else{
		int maxstate=1;
		vector<string> searchtokens;
		Tokenize(seqs[0].get_sequence(), searchtokens, " ");
		for(unsigned int j=0;j<searchtokens.size();j++){
			TrimSpaces(searchtokens[j]);
		}
		nsites = searchtokens.size();
		if(verbose)
			cout << "nsites: " << nsites << endl;
		for(unsigned int se = 0;se<seqs.size();se++){
			searchtokens = vector<string> ();
			Tokenize(seqs[se].get_sequence(), searchtokens, " 	");
			for(unsigned int j=0;j<searchtokens.size();j++){
				TrimSpaces(searchtokens[j]);
				int pos = atoi(searchtokens[j].c_str());
				if (pos > maxstate)
					maxstate = pos;
			}
		}
		nstates = maxstate+1;//TODO this can be determined by largest number +1
	}
	if(verbose)
		cout << "total number of states in dataset: " << nstates << endl;
	
	ofstream ancout;
	ofstream stnumout;
	ofstream sttimeout;
	ofstream sttnumout_any;
	
	if (ancstates.size() > 0 && outfile_anc != ""){
		ancout.open(outfile_anc.c_str(),ios::out);
		ancout << "site\ttree\tMRCA\tlnL";
		for(int i=0;i<nstates;i++){
			ancout << "\tstate_" << i+1;
		}
		ancout << endl;
	}
	if (stochnumber.size() > 0 && outfile_stochnum != ""){
		stnumout.open(outfile_stochnum.c_str(),ios::out);
		stnumout << "site\ttree\tMRCA\tlnL";
		for(int i=0;i<nstates;i++){
			for(int j=0;j<nstates;j++){
				if (i != j)
					stnumout << "\tstate_" << i+1 << "->state_" << j+1;
			}
		}
		stnumout << endl;
	}
	if (stochtime.size() > 0 && outfile_stochtime != "" ){
		sttimeout.open(outfile_stochtime.c_str(),ios::out);
		sttimeout << "site\ttree\tMRCA\tlnL";
		for(int i=0;i<nstates;i++){
			sttimeout << "\tstate_" << i+1;
		}
		sttimeout << endl;
	}
	if (stochnumber_any.size() > 0 && outfile_stochnum_any != "" ){
		sttnumout_any.open(outfile_stochnum_any.c_str(),ios::out);
		sttnumout_any << "site\ttree\tMRCA\tlnL";
		sttnumout_any << "\tanystate";
		sttnumout_any << endl;
	}
	
	for(int n = 0;n<nsites;n++){
		cout << "site: " << n+1 << endl;
		/*
		 * this converts the data and is a little long to accomodate datasets
		 * with sites that don't have all the states but the results can still
		 * be printed to the same outfile for analysis after
		 */
		//need to put the data into a wide view
		vector<Sequence> runseqs;
		int nstates_site_n;
		vector<int> existing_states(nstates,0);
		if(datawide == false){
			for(unsigned int se = 0;se<seqs.size();se++){
				vector<string> searchtokens;
				Tokenize(seqs[se].get_sequence(), searchtokens, " 	");
				for(unsigned int j=0;j<searchtokens.size();j++){
					TrimSpaces(searchtokens[j]);
				}
				string tseqs(nstates,'0');
				if(searchtokens[n]=="?"){
					for(int mse = 0;mse<nstates;mse++){
						tseqs.replace(mse,1,"1");
					}
				}else{
					int pos = atoi(searchtokens[n].c_str());
					tseqs.replace(pos,1,"1");
				}
				for(int i=0;i<nstates;i++){
					if(tseqs.at(i) =='1'){
						existing_states[i] = 1;
					}
				}
				Sequence tse =Sequence(seqs[se].get_id(),tseqs,true);
				runseqs.push_back(tse);
			}
			nstates_site_n = calculate_vector_int_sum(existing_states);
		}else{
			runseqs = seqs;
			for(unsigned int se = 0;se<seqs.size();se++){
				for(int i=0;i<nstates;i++){
					if(seqs[se].get_sequence().at(i) =='1'){
						existing_states[i] = 1;
					}
				}
			}
			nstates_site_n = calculate_vector_int_sum(existing_states);
		}
		//mapping the existing states to the full states
		int statecnt = 0;
		for(int i=nstates-1;i>=0;i--){
			if(existing_states[i] == 1){
				continue;
			}else{
				for(unsigned int se =0 ; se<runseqs.size();se++){
					runseqs[se].set_sequence(runseqs[se].get_sequence().erase(i,1));
				}
			}
		}

		if(verbose)
			cout <<"states: " << nstates_site_n << endl;
		if(verbose)
			cout << "trees: ";
		for(unsigned int i=0;i<lines.size();i++){
			if(verbose)
				cout << i << endl;
			RateModel rm(nstates_site_n);
			rm.setup_P(0.1,false);
			StateReconstructor sr(rm);
			sr.set_store_p_matrices(false);
			Tree * tree = tr.readTree(lines[i]);
			if(verbose)
				cout << "tips: "<< tree->getExternalNodeCount() << endl;

			sr.set_tree(tree);
			bool same = sr.set_tip_conditionals(runseqs);
			if (same == true){
				cout << "skipping calculation" <<endl;
				continue;
			}

			mat free_var(nstates_site_n,nstates_site_n);free_var.fill(0);
			int ct = 0;
			if(freeparams == "_one_"){
				ct = 1;
			}else if(freeparams == "_all_"){
				ct = 0;
				for(int k=0;k<nstates_site_n;k++){
					for(int j=0;j<nstates_site_n;j++){
						if(k != j){
							free_var(k,j) = ct;
							ct += 1;
						}
					}
				}
			}
			if(verbose)
				cout << free_var << endl;
			if(verbose)
				cout << ct << endl;
			rm.neg_p = false;
			optimize_sr_nlopt(&rm,&sr,&free_var,ct);
			if(verbose)
				cout << free_var << endl;
			rm.setup_Q(free_var);

			sr.set_store_p_matrices(true);
			double totlike;
			double finallike = sr.eval_likelihood();
			if(verbose)
				cout << "final_likelihood: " << finallike << endl;

			if(verbose)
				cout << "ancestral states" <<endl;
			sr.prepare_ancstate_reverse();
			for(unsigned int j=0;j<ancstates.size();j++){
				vector<double> lhoods;
				if(verbose)
					cout <<"node: " << tree->getMRCA(mrcas[ancstates[j]])->getName() << "\tmrca: " << ancstates[j] <<  endl;
				ancout << n+1 << "\t" << i+1 << "\t" << ancstates[j] << "\t" << finallike;
				lhoods = sr.calculate_ancstate_reverse(*tree->getMRCA(mrcas[ancstates[j]]));
				totlike = calculate_vector_double_sum(lhoods);
				//cout << totlike << " " << log(totlike) << endl;
				bool neg = false;
				int excount = 0;
				for(int k=0;k<nstates;k++){
					if(existing_states[k] == 1){
						if(verbose)
							cout << lhoods[excount]/totlike << " ";//"(" << lhoods[excount] << ") ";
						ancout << "\t" << lhoods[excount]/totlike;
						if (lhoods[excount]/totlike < 0)
							neg = true;
						excount += 1;
					}else{
						if(verbose)
							cout << "NA" << " ";
						ancout << "\t" << "NA";
					}
				}
				if (neg == true)
					exit(0);
				ancout <<endl;
				if(verbose)
					cout << endl;
			}
			if(verbose)
				cout << endl;

			if(verbose)
				cout << "stochastic time" << endl;
			for(unsigned int j=0;j<stochtime.size();j++){
				if(tree->getMRCA(mrcas[stochtime[j]])->isRoot() == false){
					vector<double> lhoods;
					if(verbose)
						cout <<"node: " << tree->getMRCA(mrcas[stochtime[j]])->getName() << " mrca: " << stochtime[j] <<  endl;
					sttimeout << n+1 << "\t" << i+1 << "\t" << stochtime[j]<< "\t" << finallike;
					bool neg = false;
					int excount = 0;
					for(int k=0;k<nstates;k++){
						if(existing_states[k]==1){
							sr.prepare_stochmap_reverse_all_nodes(excount,excount);
							sr.prepare_ancstate_reverse();
							vector<double> stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochtime[j]]),true);
							double tnum = calculate_vector_double_sum(stoch)/totlike;
							if(verbose)
								cout << tnum << " ";
							sttimeout << "\t" << tnum;
							if (tnum < 0)
								neg = true;
							excount += 1;
						}else{
							if(verbose)
								cout << "NA" << " ";
							sttimeout << "\t" << "NA";
						}

					}
					sttimeout << endl;
					if(verbose)
						cout << endl;
					if (neg == true)
						exit(0);
				}
			}
			if(verbose)
				cout << endl;

			if(verbose)
				cout << "stochastic number" << endl;
			for(unsigned int j=0;j<stochnumber.size();j++){
				if(tree->getMRCA(mrcas[stochnumber[j]])->isRoot() == false){
					vector<double> lhoods;
					if(verbose)
						cout <<"node: " << tree->getMRCA(mrcas[stochnumber[j]])->getName() << " mrca: " << stochnumber[j] <<  endl;
					stnumout << n+1 << "\t" << i+1 << "\t" << stochnumber[j]<< "\t" << finallike;
					bool neg = false;
					int excount = 0;
					for(int k=0;k<nstates;k++){
						if(existing_states[k]==1){
							int excount2 = 0;
							for(int l=0;l<nstates;l++){
								if(existing_states[l]==1){
									if(k==l){
										if(verbose)
											cout << " - ";
									}else{
										sr.prepare_stochmap_reverse_all_nodes(excount,excount2);
										sr.prepare_ancstate_reverse();
										vector<double> stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochnumber[j]]),false);
										double tnum = calculate_vector_double_sum(stoch)/totlike;
										if(verbose)
											cout << tnum << " ";
										stnumout << "\t" << tnum;
										if (tnum < 0)
											neg = true;
									}
									excount2 += 1;
								}else{
									if(verbose)
										cout << "NA" << " ";
									stnumout << "\t" << "NA";
								}
							}
							if(verbose)
								cout << endl;
							excount += 1;
						}else{
							for(int l=0;l<nstates;l++){
								if(k==l){
									if(verbose)
										cout << " - ";
								}else{
									if(verbose)
										cout << "NA" << " ";
									stnumout << "\t" << "NA";
								}
							}
							if(verbose)
								cout << endl;
						}
					}
					stnumout << endl;
					if(verbose)
						cout << endl;
					if (neg == true)
						exit(0);
				}
			}
			if(verbose)
				cout << endl;

			if(verbose)
				cout << "stochastic number (any)" << endl;
			if(stochnumber_any.size() > 0){
				sr.prepare_stochmap_reverse_all_nodes_all_matrices();
				sr.prepare_ancstate_reverse();
			}
			for(unsigned int j=0;j<stochnumber_any.size();j++){
				if(tree->getMRCA(mrcas[stochnumber_any[j]])->isRoot() == false){
					vector<double> lhoods;
					if(verbose)
						cout <<"node: " << tree->getMRCA(mrcas[stochnumber_any[j]])->getName() << " mrca: " << stochnumber_any[j] <<  endl;
					sttnumout_any << n+1 << "\t" << i+1 << "\t" << stochnumber_any[j]<< "\t" << finallike;

					vector<double> stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochnumber_any[j]]),false);
					double tnum = calculate_vector_double_sum(stoch)/totlike;
					//cout << calculate_vector_double_sum(stoch)<< " "<<totlike << endl;
					if(verbose)
						cout << tnum << " ";
					sttnumout_any << "\t" << tnum;
					sttnumout_any << endl;
					if(verbose)
						cout << endl;
				}
			}

			delete tree;
		}
	}
	if (ancstates.size() > 0  && outfile_anc != ""){
		ancout.close();
	}
	if (stochnumber.size() > 0){
		stnumout.close();
	}
	if (stochtime.size() > 0){
		sttimeout.close();
	}
	
	return EXIT_SUCCESS;
}
