#ifndef PX_NJ_H
#define PX_NJ_H

#include <map>
#include <vector>
#include <string>
#include <iostream>

class NJOI {
private:
    std::map<std::string, std::string> sequences_;
    std::map<std::string, std::string>::iterator iter_;
    std::vector<std::string> names_;
    //std::string fasta;
    std::vector< std::vector<double> > Matrix_;
    void set_name_key ();

    // additions:
    int num_taxa_;
    int num_char_;
    int nthreads_; // not implemented yet
    std::map<int, std::string> name_key_;
    std::string newick_string_; // temporary
    
    void CalcQ (const int& NumbOfSequences, std::vector< std::vector<double> >& OriginalMatrix, 
        std::vector< std::vector<double> >& ConvertedMatrix, std::vector< std::vector<double> >& LengthMatrix);
    void FetchLengths (const int& NumbOfSequences, const std::vector< std::vector<double> >& NewMatrix,
        std::vector< std::vector<double> >& LengthMatrix, const int& mini1, const int& mini2,
        double& brlength1, double& brlength2);
    void Tree_Update (std::string& newname, std::vector<std::string>& names,
        int& NumbOfSequences, std::vector< std::vector<double> >& NewMatrix, int& mini1, int& mini2,
        double& brlength1, double& brlength2);
    void Choose_Smallest (int& NumbOfSequences, const std::vector< std::vector<double> >& Matrix,
        int& mini1, int& mini2);
    
public:
    NJOI (std::istream* pios, int& threads);
    std::map<std::string, std::string> FastaToOneLine (std::string& fasta);
    std::vector< std::vector<double> > BuildMatrix (std::map<std::string, std::string>& sequences);
    void TREEMAKE (std::vector<std::string>&, std::map<int, std::string>&,
        std::vector< std::vector<double> >&);
    std::string get_newick ();
};

#endif /* PX_NJ_H */
