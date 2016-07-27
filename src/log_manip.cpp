#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <numeric>

using namespace std;

#include "log_manip.h"
#include "utils.h"

// this doesn't quite work with an optional arg for input_files
LogManipulator::LogManipulator(string const& logtype, vector <string> const& input_files,
    istream* pios, ostream* poos) {
    if (input_files.size() > 0) {
        files_ = input_files;
        num_files_ = input_files.size();
    } else {
        pios_ = pios;
    }
    logtype_ = logtype;
    poos_ = poos;
}

LogManipulator::LogManipulator(string const& logtype, vector <string> const& input_files,
    ostream* poos) {
    files_ = input_files;
    num_files_ = input_files.size();
    logtype_ = logtype;
    poos_ = poos;
}

void LogManipulator::sample(int const& burnin, int const& nthin, int const& nrandom,
    int const& seed) {
    burnin_ = burnin;
    nthin_= nthin;
    nrandom_ = nrandom;
    seed_ = seed;
    if (logtype_ == "parm") {
        sample_parameters ();
    } else {
        sample_trees ();
    }
}

void LogManipulator::count () {
    if (logtype_ == "parm") {
        count_parameter_samples ();
    } else {
        count_tree_samples ();
    }
}

void LogManipulator::count_parameter_samples () {
    num_cols_ = 0;
    if (!files_.empty()) {
        for (int i=0; i < num_files_; i++) {
            string curfile = files_[i];
            infilestr_.open(curfile.c_str());
            string line;
            bool first_line = true;
            int num_samps = 0;
            while (getline(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                } else if (first_line) {
                    vector <string> header = tokenize(line);
                    int curpars = header.size();
                    if (i == 0) { // first header
                        num_cols_ = curpars;
                        parm_columns_ = header;
                    } else {
                        // check that we've still got the same number of parameters i.e. files match
                        if (curpars != num_cols_) {
                            cout << "Error: number of parameters in file " << (i + 1)
                                << "(" << curpars << ") does not match that from first file ("
                                << num_cols_ << "). Exiting." << endl;
                            exit(0);
                        } else if (header != parm_columns_) {
                            // check that headers are identical
                            cout << "Error: header for file " << (i + 1)
                                << "does not match that from first file. Exiting." << endl;
                            exit(0);
                        }
                    }
                    
                    first_line = false;
                    continue;
                } else {
                    num_samps++;
                    continue;
                }
            }
            indiv_totals_.push_back(num_samps);
            infilestr_.close();
        }
        ntotal_samples_ = accumulate(indiv_totals_.begin(), indiv_totals_.end(), 0);
        (*poos_) << "Counted " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
            << " variables across " << num_files_ << " files." << endl;
    } else {
        
        // stream stuff will go here (maybe)
        
    }
    
}

void LogManipulator::count_tree_samples () {
    
}

void LogManipulator::sample_parameters () {
    if (!files_.empty()) {
        for (int i=0; i < num_files_; i++) {
            string curfile = files_[i];
            infilestr_.open(curfile.c_str());
            string line;
            bool first_line = true;
            int num_samps = 0;
            while (getline(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                } else if (first_line) {
                    vector <string> header = tokenize(line);
                    int curpars = header.size();
                    if (i == 0) { // first header
                        num_cols_ = curpars;
                        parm_columns_ = header;
                        for (int j=0; j < num_cols_; j++) {
                            (*poos_) << parm_columns_[j];
                            if (j < (num_cols_ - 1)) {
                                (*poos_) << "\t";
                            }
                        }
                        (*poos_) << endl;
                    } else {
                        // check that we've still got the same number of parameters i.e. files match
                        if (curpars != num_cols_) {
                            cout << "Error: number of parameters in file " << (i + 1)
                                << "(" << curpars << ") does not match that from first file ("
                                << num_cols_ << "). Exiting." << endl;
                            exit(0);
                        } else if (header != parm_columns_) {
                            // check that headers are identical
                            cout << "Error: header for file " << (i + 1)
                                << "does not match that from first file. Exiting." << endl;
                            exit(0);
                        }
                    }
                    first_line = false;
                    continue;
                } else {
                    if ((num_samps - burnin_) > 0 && (num_samps - burnin_) < nthin_) {
                        num_samps++;
                        (*poos_) << line << endl;
                    }
                    continue;
                }
            }
            indiv_totals_.push_back(num_samps);
            infilestr_.close();
        }
        ntotal_samples_ = accumulate(indiv_totals_.begin(), indiv_totals_.end(), 0);
        (*poos_) << "Counted " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
            << " variables across " << num_files_ << " files." << endl;
    }
}

void LogManipulator::sample_trees () {
    
}
