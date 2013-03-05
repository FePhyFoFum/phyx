#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iostream>
#include <fstream>

using namespace std;

#include "seq_models.h"
#include "sequence.h"
#include "utils.h"

void read_scoring_matrix(char * filename, map<char, map<char,int> > & sc_mat){
    ifstream fstr(filename);
    sc_mat.clear();
    vector<char> order;
    bool first = true;
    string line;
    while(getline(fstr,line)){
        if (line[0] == '#'){
            continue;
        }else{
            vector<string> tokens;
            string del(" \t");
            tokenize(line,tokens,del);
            for(int i=0;i<tokens.size();i++){
                trim_spaces(tokens[i]);
            }
            if (first == true){
                first = false;
                for(int i=0;i<tokens.size();i++){
                    order.push_back(tokens[i][0]);
                }
                for (int j=0;j<order.size();j++){
                    sc_mat[order[j]] = map<char,int>();
                }
                continue;
            }
            for (int j=0;j<order.size();j++){
                sc_mat[tokens[0][0]][order[j]] = atoi(tokens[j+1].c_str()); //#changed from int to float
            }
        }
    }
    fstr.close();
}


