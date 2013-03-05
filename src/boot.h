#ifndef SEQ_SAMP_H_
#define SEQ_SAMP_H_

class SequenceSampler {
private:
    float jkfract;
    bool jackknife;
    vector <int> samplesites;
    
    void get_bootstrap_sites (int const&);
    void get_jackknife_sites (int const&);

public:
    SequenceSampler (int const&, float const&);
    vector <int> get_sampled_sites ();
    void sample_sites (int const&);
    string get_resampled_seq (string const& origseq);
    //~SequenceResampler();
};

#endif /* SEQ_SAMP_H_ */
