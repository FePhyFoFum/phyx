
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iostream>

using namespace std;

#include "seq_utils.h"
#include "sequence.h"

/**
 * TODO: need to make a alphabet guesser
 */

/**
 * procedure to get the character of a nucleotide
 * from the set of positions
 */
char get_dna_from_pos(set<int> ins){
    if(ins.count(0) == 1){
        if(ins.count(1) == 1){
            if(ins.count(2) == 1){
                if(ins.count(3) == 1){
                    return 'N';
                }
                return 'V';
            }
            if(ins.count(3) == 1){
                return 'H';
            }
            return 'M';
        }
        if(ins.count(2) == 1){
            return 'R';
        }
        if(ins.count(3) == 1){
            return 'W';
        }
        return 'A';
    }
    if(ins.count(1)==1){
        if(ins.count(2) == 1){
            if(ins.count(3) == 1){
                return 'B';
            }
            return 'S';
        }
        if(ins.count(3) == 1){
            return 'Y';
        }
        return 'C';
    }
    if(ins.count(2)==1){
        if(ins.count(3) == 1){
            return 'K';
        }
        return 'G';
    }
    if(ins.count(3)==1){
        return 'T';
    }
}

set<int> get_dna_pos(char inc){
    set<int> ret;
    inc = toupper(inc);
    if(inc == 'A'){
        ret.insert(0);
    }else if(inc == 'C'){
        ret.insert(1);
    }else if(inc == 'G'){
        ret.insert(2);
    }else if(inc == 'T'){
        ret.insert(3);
    }else if(inc == '-' || inc == 'N'){
        ret.insert(0);ret.insert(1);ret.insert(2);ret.insert(3);
    }else if(inc == 'Y'){
        ret.insert(1);ret.insert(3);
    }else if(inc == 'R'){
        ret.insert(0);ret.insert(2);
    }else if(inc == 'W'){
        ret.insert(0);ret.insert(3);
    }else if(inc == 'M'){
        ret.insert(0);ret.insert(1);
    }else if(inc == 'B'){
        ret.insert(1);ret.insert(2);ret.insert(3);
    }else if(inc == 'V'){
        ret.insert(0);ret.insert(1);ret.insert(2);
    }else if(inc == 'S'){
        ret.insert(1);ret.insert(2);
    }else if(inc == 'K'){
        ret.insert(2);ret.insert(3);
    }else if(inc == 'H'){
        ret.insert(0);ret.insert(1);ret.insert(3);
    }
    return ret;
}

/**
 * 
 * int alpha: the alphabet with 0=dna, 1=aa
 */
string consensus_seq(vector<Sequence> & seqs, int alpha){
    int seqlength = seqs[0].get_sequence().length();
    for(int i=0;i<seqs.size();i++){
        assert(seqs[i].get_sequence().length() == seqlength);
    }
    string retstring;
    for(int i=0;i<seqlength;i++){
        set<int> fullset;
        for(int j=0;j<seqs.size();j++){
            set<int> tset = get_dna_pos(seqs[j].get_sequence()[i]);
            fullset.insert(tset.begin(),tset.end());
        }
        retstring += get_dna_from_pos(fullset);
    }
    return retstring;
}
