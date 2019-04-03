#include <vector>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <algorithm>

using namespace std;

#include "sequence.h"
#include "pairwise_alignment.h"

/**
 * Needleman-Wunsch
 * returning aln1, aln2 strings and score
 *
 * scoringmatrix should come from read_score_matrix
 */
double nw(Sequence & iseq1, Sequence & iseq2, map<char, map<char,int> > & scoringmatrix,
    double gap_penalty, string & aln1, string & aln2) {
    
    string seq1 = iseq1.seq_to_upper();
    string seq2 = iseq2.seq_to_upper();
    
    vector<vector<double> > F;
    for (unsigned int i=0; i < seq2.length()+1; i++) {
        vector<double> b;
        for (unsigned int j=0; j < seq1.length()+1; j++) {
            b.push_back(0);
        }
        F.push_back(b);
    }
    double d = gap_penalty; // just easier to read
    for (unsigned int i=0; i < seq1.length()+1; i++) {
        F[0][i] = d*i;
    }
    for (unsigned int j=0; j < seq2.length()+1; j++) {
        F[j][0] = d*j;
    }
    for (unsigned int j=1; j < seq2.length()+1; j++) {
        for (unsigned int i=1; i < seq1.length()+1; i++) {
            double mat = F[j-1][i-1] + scoringmatrix[seq1[i-1]][seq2[j-1]];
            double del = F[j-1][i] + d;
            double ins = F[j][i-1] + d;
            double v = mat;
            if (del > v) {
                v = del;
            }
            if (ins > v) {
                v = ins;
            }
            F[j][i] = v;
        }
    }
    aln1 = "";
    aln2 = "";
    int i = seq1.length();
    int j = seq2.length();
    while (i > 0 || j > 0) {
        double score = F[j][i];
        double scorediag = -999;
        if (j > 0 && i > 0) {
            scorediag = F[j-1][i-1];
        }
        double scoreup = -999;
        if (j > 0) {
            scoreup = F[j-1][i];
        }
                double scoreleft = -999;
        if (i > 0) {
            scoreleft = F[j][i-1];
        }
        if (i > 0 && j > 0 && score == scorediag + scoringmatrix[seq1[i-1]][seq2[j-1]]) {
            aln1.append(1,seq1[i-1]);
            aln2.append(1,seq2[j-1]);
            i = i - 1;
            j = j - 1;
        } else if ( i > 0 && score == scoreleft + d) {
            aln1.append(1,seq1[i-1]);
            aln2.append("-");
            i = i - 1;
        } else if (j > 0 && score == scoreup + d) {
            aln1.append("-");
            aln2.append(1,seq2[j-1]);
            j = j - 1;
        }
    }
    reverse(aln1.begin(),aln1.end());
    reverse(aln2.begin(),aln2.end());
    double score=0;
    for (unsigned int i=0; i < aln1.length(); i++) {
        if (aln1[i]!='-' && aln2[i]!='-') {
            score=score+scoringmatrix[aln1[i]][aln2[i]];
        } else {
            score=score+gap_penalty;
        }
    }
    return score;
}

/**
 * Smith-Waterman
 * returning aln1, aln2 strings and score
 *
 * scoringmatrix should come from read_score_matrix
 */
double sw(Sequence & iseq1, Sequence & iseq2, map<char, map<char,int> > & scoringmatrix,
    double gap_penalty, string & aln1, string & aln2) {
    
    string seq1 = iseq1.seq_to_upper();
    string seq2 = iseq2.seq_to_upper();
    
    vector<vector<double> > F;
    for (unsigned int i=0; i < seq2.length()+1; i++) {
        vector<double> b;
        for (unsigned int j=0; j < seq1.length()+1; j++) {
            b.push_back(0);
        }
        F.push_back(b);
    }
    double d = gap_penalty; // just easier to read
    for (unsigned int i=0; i < seq1.length()+1; i++) {
        F[0][i] = d*i;
    }
    for (unsigned int j=0; j < seq2.length()+1; j++) {
        F[j][0] = d*j;
    }
    int besti = 0;
    int bestj = 0;
    double bestsc = 0;
    for (unsigned int j=1; j < seq2.length()+1; j++) {
        for (unsigned int i=1; i < seq1.length()+1; i++) {
            double mat = F[j-1][i-1] + scoringmatrix[seq1[i-1]][seq2[j-1]];
            double del = F[j-1][i] + d;
            double ins = F[j][i-1] + d;
            double v = mat;
            if (del > v) {
                v = del;
            }
            if (ins > v) {
                v = ins;
            }
            F[j][i] = v;
            if (v > bestsc) {
                bestsc = v;
                besti = i;
                bestj = j;
            }
        }
    }    
    aln1 = "";
    aln2 = "";
    int i = besti;
    int j = bestj;
    while (i > 0 && j > 0) {
        double score = F[j][i];
        if (score == 0) {
                break;
        }
        double scorediag = -999;
        if ( j > 0 && i > 0) {
            scorediag = F[j-1][i-1];
        }
        double scoreup = -999;
        if (j > 0) {
            scoreup = F[j-1][i];
        }
        double scoreleft = -999;
        if (i > 0) {
            scoreleft = F[j][i-1];
        }
        if (i > 0 && j > 0 && score == scorediag + scoringmatrix[seq1[i-1]][seq2[j-1]]) {
            aln1.append(1,seq1[i-1]);
            aln2.append(1,seq2[j-1]);
            i = i - 1;
            j = j - 1;
        } else if (i > 0 && score == scoreleft + d) {
            aln1.append(1,seq1[i-1]);
            aln2.append("-");
            i = i - 1;
        } else if (j > 0 && score == scoreup + d) {
            aln1.append("-");
            aln2.append(1,seq2[j-1]);
            j = j - 1;
        }
    }
    reverse(aln1.begin(),aln1.end());
    reverse(aln2.begin(),aln2.end());
    double score=0;
    for (unsigned int i=0; i < aln1.length(); i++) {
        if (aln1[i]!='-' && aln2[i]!='-') {
            score=score+scoringmatrix[aln1[i]][aln2[i]];
        } else {
            score=score+gap_penalty;
        }
    }
    return score;
}
