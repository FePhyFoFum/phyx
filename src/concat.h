#ifndef _CONCAT_H_
#define _CONCAT_H_

#include <vector>

using namespace std;

#include "sequence.h"

class SequenceConcatenater {
private:
    vector <Sequence> seqs_;
    int num_partitions_;
    int num_char_;
    int num_taxa_;
    int ft_;
    bool toupcase_;
    bool interleave_;
    string filename_;
    vector <int> partition_sizes_;
    
    void read_sequences (string & seqf);
    void delete_sequence (SequenceConcatenater & newSeqs, int const& index);

public:
    SequenceConcatenater ();
    SequenceConcatenater (string & seqf, bool & toupcase);
    void concatenate (SequenceConcatenater & newSeqs);
    int get_sequence_length ();
    int get_num_taxa ();
    Sequence get_sequence (int const & index);
    vector <int> get_partition_sizes ();
    void write_partition_information (vector <string> const& inputFiles,
        string & partfile);
    //~SequenceConcatenater();
};

#endif /* _CONCAT_H_ */
