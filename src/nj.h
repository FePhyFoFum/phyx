#ifndef _NJ_H_
#define _NJ_H_

#include <map>
#include <vector>
#include <string>

using namespace std;

class NJOI {
private:

    map <string, string> sequences_;
    map <string, string>::iterator iter_;
    vector<string> names_;
    //string fasta;
    vector< vector<double> > Matrix;
    void set_name_key ();

    // additions:
    int ntax_;
    int nchar_;
    int nthreads_; // not implemented yet
    map<int, string> name_key_;
    string newick_string_; // temporary
    
    void CalcQ(int const& NumbOfSequences, vector< vector<double> >& OriginalMatrix, 
        vector< vector<double> >& ConvertedMatrix, vector< vector<double> >& LengthMatrix);
    void FetchLengths(int const& NumbOfSequences, vector< vector<double> > const& NewMatrix,
    vector< vector<double> >& LengthMatrix, int const& mini1, int const& mini2,
        double & brlength1, double & brlength2);
    void Tree_Update(string& newname, vector<string>& names, map<int, string>& NumbKeys,
        int& NumbOfSequences, vector< vector<double> >& NewMatrix, int& mini1, int& mini2,
        double& brlength1, double& brlength2);
    void Choose_Smallest(int& NumbOfSequences, vector< vector<double> > const& Matrix,
        int & mini1, int & mini2);
    
public:
    NJOI (istream* pios, int & threads);
    map<string, string> FastaToOneLine(string& fasta);
    vector< vector<double> > BuildMatrix(map<string, string>& sequences);
    void TREEMAKE(vector<string>&, map <int, string>&, vector< vector<double> >&);
    string get_newick ();
};

#endif /* _NJ_H_ */
