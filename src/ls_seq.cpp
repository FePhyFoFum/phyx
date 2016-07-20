#include <string>
#include <map>
#include <iomanip>
#include <iostream>

using namespace std;

#include "ls_seq.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

void Stats::STAT_Getter(string& seq, bool& prot){

    Total.clear();
    for (int i = 0; i < Molecule.length(); i++){
        Total[Molecule[i]] = 0.0;
    }
    for (int i = 0; i < seq.length(); i++){
        seq[i] = toupper(seq[i]);
        //Ensure theres no weird J or whatever characters
        if (Total.find(seq[i]) == Total.end()) {
            if (prot == true){
                Total['X']++;
            }else{
                Total['N']++;
            }
        }else{
            Total[seq[i]]++;
        }
    }
}

// *** this should not print to cout, but instead to pios ***
void Stats::Printer(bool& prot){
    
        const char separator = ' ';
        const int nameWidth = 10;
        double divide = 0.0;
        if (prot == true){
            Mol = "Prot ";
        }else{
            Mol = "Nucl ";
        }
        if (finished == true){
            cout << "General Stats For All Sequences" << endl;
            cout << "File Type: " << type << endl;
            cout << "Number of Sequences: " << seqcount << endl;
            cout << "Total Length of All Combined: " << Concatenated.length() << endl;
            divide = Concatenated.length();
        }else{
            cout << "General Stats For " << name << endl;
            cout << "Total Length: " << temp_seq.length() << endl;    
            divide = temp_seq.length();
        }
        cout << "--------" << Mol << "TABLE---------" << endl;
        cout << Mol << "\tTotal\tPercent" << endl;
        for (int i = 0; i < Molecule.length(); i++){
            cout << left << setw(nameWidth) << setfill(separator) << Molecule[i] << Total[Molecule[i]] << "\t" << ((Total[Molecule[i]] / divide)*100.0) << endl;
        }
        if (prot == false){
        cout << left << setw(nameWidth) << setfill(separator) << "G+C" << (Total['G'] + Total['C'])<< "\t" << (((Total['G'] + Total['C']) / divide)*100.0) << endl;
    }        
    cout << "--------" << Mol << "TABLE---------" << endl;

}

Stats::Stats (istream* pios, bool& all, bool& prot) {

    //Concatenated will be used for all stats
    finished = false;
    seqcount = 0;
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    if (prot == true){
        Molecule = "ACDEFGHIKLMNPQRSTVWXY*";
    }else{
        Molecule = "ACGTN-";
    }
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        
        seqcount++;
        Concatenated += seq.get_sequence();
        temp_seq = seq.get_sequence();
        name = seq.get_id();
        if (all == true){
            STAT_Getter(temp_seq, prot);
            Printer(prot);
        }
        if (ft == 1){
            type = "Phylip";
        }
        if (ft == 0){
            type = "Nexus";
        }
    
    }
    if (ft == 2) {
        seqcount++;
        Concatenated += seq.get_sequence();
        temp_seq = seq.get_sequence();
        name = seq.get_id();
        type = "Fasta";
        if (all == true){
            STAT_Getter(temp_seq, prot);
            Printer(prot);
        }
    }
    finished = true;
    STAT_Getter(Concatenated, prot);
    Printer(prot);
}
