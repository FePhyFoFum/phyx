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
// so doesn't currently read from stream (because need to specify parameter vs.tree)
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
    ostream* poos, bool const& verbose) {
    files_ = input_files;
    num_files_ = input_files.size();
    logtype_ = logtype;
    poos_ = poos;
    verbose_ = verbose;
}

void LogManipulator::sample(int const& burnin, int const& nthin, int const& nrandom,
    int const& seed) {
    burnin_ = burnin;
    nthin_= nthin;
    nrandom_ = nrandom;
    seed_ = seed;
    if (logtype_ == "parameter") {
        sample_parameters ();
    } else {
        sample_trees ();
    }
}

void LogManipulator::count () {
    if (logtype_ == "parameter") {
        count_parameter_samples ();
    } else {
        count_tree_samples ();
    }
}

// TODO: should counting allow burnin/thinning?
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
            indiv_sample_totals_.push_back(num_samps);
            infilestr_.close();
        }
        //ntotal_samples_ = accumulate(indiv_totals_.begin(), indiv_totals_.end(), 0);
        //(*poos_) << "Counted " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
        //    << " variables across " << num_files_ << " files." << endl;
    } else {
        
        // stream stuff will go here (maybe)
        
    }
    ntotal_samples_ = accumulate(indiv_sample_totals_.begin(), indiv_sample_totals_.end(), 0);
}

void LogManipulator::count_tree_samples () {
    if (!files_.empty()) {
        for (int i=0; i < num_files_; i++) {
            string curfile = files_[i];
            infilestr_.open(curfile.c_str());
            string line;
            int num_samps = 0;
            while (getline(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                } else {
                    vector <string> tokens = tokenize(line);
                    string first = tokens[0];
                    std::transform(first.begin(), first.end(), first.begin(), ::tolower);
                    if (first == "tree") {
                        num_samps++;
                    }
                }
            }
            indiv_sample_totals_.push_back(num_samps);
            infilestr_.close();
        }
        //ntotal_samples_ = accumulate(indiv_totals_.begin(), indiv_totals_.end(), 0);
        //(*poos_) << "Counted " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
        //    << " variables across " << num_files_ << " files." << endl;
    } else {
        
        // stream stuff will go here (maybe)
        
    }
    ntotal_samples_ = accumulate(indiv_sample_totals_.begin(), indiv_sample_totals_.end(), 0);
}

void LogManipulator::get_sample_counts () {
    if (logtype_ == "parameter") {
        for (int i = 0; i < num_files_; i++) {
            (*poos_) << files_[i] << ": " << indiv_sample_totals_[i] << " samples of "
                << (num_cols_ - 1) << " variables." << endl;
        }
        if (num_files_ > 1) {
            (*poos_) << "Counted " << ntotal_samples_ << " total samples of "
                << (num_cols_ - 1) << " variables across " << num_files_ << " files." << endl;
        }
    } else {
        for (int i = 0; i < num_files_; i++) {
            (*poos_) << files_[i] << ": " << indiv_sample_totals_[i] << " trees." << endl;
        }
        if (num_files_ > 1) {
            (*poos_) << "Counted " << ntotal_samples_ << " total tree samples across "
                << num_files_ << " files." << endl;
        }
    }
}

void LogManipulator::get_column_names () {
    if (!files_.empty()) {
        if (num_files_ == 1) {
            string curfile = files_[0];
            infilestr_.open(curfile.c_str());
            string line;
            bool first_line = true;
            while (getline(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                } else if (first_line) {
                    vector <string> header = tokenize(line);
                    int curpars = header.size();
                    num_cols_ = curpars;
                    parm_columns_ = header;
                    first_line = false;
                    break;
                }
            }
            for (int i = 0; i < num_cols_; i++) {
                (*poos_) << i+1 << ". " << parm_columns_[i] << endl;
            }
        } else {
        // multiple files. make sure number/name/order of columns is identical
            
        }
    }
}

// delete variables (columns) from a large parameter log
// assumes a single log file
// assume that burnin etc. is done separately?
void LogManipulator::delete_columns (vector <int> const& col_ids) {
    if (!files_.empty()) {
        ntotal_samples_ = 0;
        string curfile = files_[0];
        infilestr_.open(curfile.c_str());
        string line;
        bool first_line = true;
        int sample_counter = 0;
        vector <int> cols_to_retain;
        
        while (getline(infilestr_, line)) {
            if (line.empty() || check_comment_line(line)) {
                continue;
            } else if (first_line) {
                vector <string> header = tokenize(line);
                int curpars = header.size();
                num_cols_ = curpars;
                parm_columns_ = header;
                
                cols_to_retain.resize(num_cols_);
                
                // check that right end of col_ids is valid (0 already checked upstream)
                // vector has been sorted
                if (col_ids.back() > num_cols_) {
                    cout << "Warning: column numbers are 1-indexed. Exiting." << endl;
                    exit (0);
                }
                
                iota(cols_to_retain.begin(), cols_to_retain.end(), 0);
                
                // remove unwanted column indices (reverse order)
                for (vector<int>::const_reverse_iterator i = col_ids.rbegin(); i < col_ids.rend(); i++) {
                    // subtract 1 because input is 1-indexed
                    cols_to_retain.erase(cols_to_retain.begin()+(*i)-1);
                }
                
                num_cols_retain_ = cols_to_retain.size();
                for (int i=0; i < num_cols_retain_; i++) {
                    (*poos_) << header[cols_to_retain[i]];
                    if (i < (num_cols_retain_ - 1)) {
                        (*poos_) << "\t";
                    }
                }
                (*poos_) << endl;
                first_line = false;
                continue;
            } else {
                vector <string> samp = tokenize(line);
                for (int i=0; i < num_cols_retain_; i++) {
                    (*poos_) << samp[cols_to_retain[i]];
                    if (i < (num_cols_retain_ - 1)) {
                        (*poos_) << "\t";
                    }
                }
                (*poos_) << endl;
                sample_counter++;
            }
        }
            
        if (verbose_) {
            cout << "Retained " << sample_counter << " samples for "
                << num_cols_retain_ << " variables." << endl;
        }
    }
}

// complement of above: retain columns passed in
void LogManipulator::retain_columns (vector <int> const& col_ids) {
    if (!files_.empty()) {
        ntotal_samples_ = 0;
        string curfile = files_[0];
        infilestr_.open(curfile.c_str());
        string line;
        bool first_line = true;
        int sample_counter = 0;
        vector <int> cols_to_retain;
        for (int i = 0; i < (int)col_ids.size(); i++) {
            // subtract 1 because input is 1-indexed
            cols_to_retain.push_back(col_ids[i] - 1);
        }
        num_cols_retain_ = cols_to_retain.size();
        
        while (getline(infilestr_, line)) {
            if (line.empty() || check_comment_line(line)) {
                continue;
            } else if (first_line) {
                vector <string> header = tokenize(line);
                int curpars = header.size();
                num_cols_ = curpars;
                
                // check that right end of col_ids is valid (0 already checked upstream)
                // vector has been sorted
                if (col_ids.back() > num_cols_) {
                    cout << "Warning: column numbers are 1-indexed. Exiting." << endl;
                    exit (0);
                }
                
                parm_columns_ = header;
                for (int i=0; i < num_cols_retain_; i++) {
                    (*poos_) << header[cols_to_retain[i]];
                    if (i < (num_cols_retain_ - 1)) {
                        (*poos_) << "\t";
                    }
                }
                (*poos_) << endl;
                first_line = false;
                continue;
            } else {
                vector <string> samp = tokenize(line);
                for (int i=0; i < num_cols_retain_; i++) {
                    (*poos_) << samp[cols_to_retain[i]];
                    if (i < (num_cols_retain_ - 1)) {
                        (*poos_) << "\t";
                    }
                }
                (*poos_) << endl;
                sample_counter++;
            }
        }
            
        if (verbose_) {
            cout << "Retained " << sample_counter << " samples for "
                << num_cols_retain_ << " variables." << endl;
        }
    }
}

// not yet used?
void LogManipulator::sample_parameters () {
    if (!files_.empty()) {
        ntotal_samples_ = 0;
        for (int i=0; i < num_files_; i++) {
            string curfile = files_[i];
            infilestr_.open(curfile.c_str());
            string line;
            bool first_line = true;
            int par_counter = 0; // this is the raw number of parameter lines in a file
            int sample_counter = 0;
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
                    if ((par_counter - burnin_) > 0 && (par_counter - burnin_) < nthin_) {
                        // skip because does not match sampling parameters
                        par_counter++;
                        continue;
                    } else if ((par_counter - burnin_) == 0) {
                        // keep first post-burnin sample from a file
                        par_counter++;
                        write_reformatted_sample(line, ntotal_samples_);
                        sample_counter++;
                        ntotal_samples_++;
                        continue;
                    } else if ((par_counter - burnin_) > 0 && (par_counter - burnin_) % nthin_ == 0) {
                        par_counter++;
                        write_reformatted_sample(line, ntotal_samples_);
                        sample_counter++;
                        ntotal_samples_++;
                        continue;
                    } else {
                        // skip because have not yet exceeded burnin
                        par_counter++;
                    }
                }
            }
            indiv_raw_counts_.push_back(par_counter);
            indiv_sample_totals_.push_back(sample_counter);
            infilestr_.close();
        }
        if (verbose_) {
            for (int i = 0; i < num_files_; i++) {
                cout << files_[i] << ": " << indiv_sample_totals_[i]
                    << " samples retained (from original " << indiv_raw_counts_[i]
                    << " samples) for " << (num_cols_ - 1) << " variables." << endl;
            }
            cout << "Retained " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
                << " variables across " << num_files_ << " input files." << endl;
        }
    }
}

void LogManipulator::sample_trees () {
    if (!files_.empty()) {
        ntotal_samples_ = 0;
        for (int i=0; i < num_files_; i++) {
            string curfile = files_[i];
            infilestr_.open(curfile.c_str());
            string line;
            int tree_counter = 0; // this is the raw number of tree lines in a file
            int sample_counter = 0;
            bool trees_encountered = false;
            while (getline(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    if (i == 0) {
                        // keep comment information from top of first file
                        (*poos_) << line << endl;
                    }
                    continue;
                } else {
                    vector <string> tokens = tokenize(line);
                    string first = tokens[0];
                    std::transform(first.begin(), first.end(), first.begin(), ::tolower);
                    if (first == "tree") {
                        if (tree_counter == 0) {
                            trees_encountered = true;
                        }
                        if (ntotal_samples_ == 0) {
                            // grab tree naming scheme
                            get_tree_name_prefix(line);
                        }
                        if ((tree_counter - burnin_) > 0 && (tree_counter - burnin_) < nthin_) {
                            // skip because does not match sampling parameters
                            tree_counter++;
                            continue;
                        } else if ((tree_counter - burnin_) == 0) {
                            // keep first post-burnin sample from a file
                            tree_counter++;
                            write_reformatted_sample(line, ntotal_samples_);
                            sample_counter++;
                            ntotal_samples_++;
                            continue;
                        } else if ((tree_counter - burnin_) > 0 && (tree_counter - burnin_) % nthin_ == 0) {
                            tree_counter++;
                            write_reformatted_sample(line, ntotal_samples_);
                            sample_counter++;
                            ntotal_samples_++;
                            continue;
                        } else {
                            // skip because have not yet exceeded burnin
                            tree_counter++;
                        }
                    } else { // not a tree line. only care about first file here. don't want anything below trees
                        // keep header from first file
                        // likely includes translation table
                        if (i == 0 && !trees_encountered) {
                            (*poos_) << line << endl;
                        }
                    }
                }
            }
            indiv_raw_counts_.push_back(tree_counter);
            indiv_sample_totals_.push_back(sample_counter);
            infilestr_.close();
        }
        (*poos_) << "End;" << endl;
        
        if (verbose_) {
            for (int i = 0; i < num_files_; i++) {
                cout << files_[i] << ": " << indiv_sample_totals_[i]
                    << " tree samples retained (from original " << indiv_raw_counts_[i]
                    << " samples)." << endl;
            }
            (*poos_) << "Retained " << ntotal_samples_ << " total tree samples across "
                << num_files_ << " input files." << endl;
        }
    }
}

// gen.NNN from Mrbayes, STATE_NNN from BEAST
// try to be general here, though: take whatever precedes the sample number
void LogManipulator::get_tree_name_prefix (string & sample) {
    vector <string> terp = tokenize(sample);
    string tree_name = terp[1];
    std::size_t found = tree_name.find_first_of("0123456789");
    tree_name.replace(tree_name.begin(), tree_name.end(), tree_name.begin(), tree_name.begin()+found);
    tree_name_prefix_ = tree_name;
}

void LogManipulator::write_reformatted_sample (string & sample, int & sample_num) {
    vector <string> terp = tokenize(sample);
    if (logtype_ == "parameter") {
        (*poos_) << sample_num;
        for (int i = 1; i < num_cols_; i++) {
            (*poos_) << "\t" << terp[i];
        }
        (*poos_) << endl;
    } else {
        // format should be: tree tree_name lots_of_other_optional_things
        (*poos_) << "tree " << tree_name_prefix_ << sample_num;
        for (unsigned int i = 2; i < terp.size(); i++) {
            (*poos_) << " " << terp[i];
        }
        (*poos_) << endl;
    }
}
