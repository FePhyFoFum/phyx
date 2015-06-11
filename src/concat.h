#ifndef _CONCAT_H_
#define _CONCAT_H_

class SequenceConcatenater {
private:
    vector <Sequence> seqs;
    int numPartitions;
    int numChar;
    int numTaxa;
    int ft;
    string filename;
    vector <int> partitionSizes;
    
    void read_sequences (string & seqf);
    void delete_sequence (SequenceConcatenater & newSeqs, int const& index);

public:
    SequenceConcatenater ();
    SequenceConcatenater (string & seqf);
    void concatenate (SequenceConcatenater & newSeqs);
    int get_sequence_length ();
    int get_num_taxa ();
    Sequence get_sequence (int const & index);
    vector <int> get_partition_sizes ();
    void write_partition_information (string & partfile);
    //~SequenceConcatenater();
};

#endif /* _CONCAT_H_ */
