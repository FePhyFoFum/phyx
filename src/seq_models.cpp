
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iostream>

using namespace std;

#include "seq_models.h"
#include "sequence.h"

void read_scoring_matrix(filename, map<char, map<char,int> > & sc_mat){
	infile = open(filename,"r")
	sc_mat = {}
	order = []
	first = True
	for i in infile:
		if i[0] == "#":
			continue
		else:
			spls = i.strip().split()
			if first == True:
				first = False
				order = spls
				for j in order:
					sc_mat[j] = {}
				continue
			for j in range(len(order)):
				sc_mat[spls[0]][order[j]] = float(spls[j+1]) #changed from int to float
	infile.close()
}
