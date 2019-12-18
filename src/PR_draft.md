# Changes to Phyx

### Repo-wide
- the biggest thing you'll notice is the removal of `using namespace std;` throughout (including header files :/). this is to avoid potential namespace conflicts. instead, things have been updated to use explicit namespace qualifications. so for example, use `std::vector<std::string>` and `std::cout` instead of `vector<string>` and `cout`. this involves more typing, but is safer. and honestly you get used to it. 
    - if this is too awkward, use an explicit (and limited) using-declaration like `using std::vector;` and `using std::cout;` to allow easy typing of `vector<int>` and `cout`, but limit the use of these and their scope as much as possible
- more extensive alignment reading capabilities, and better error-checking. this includes annoying variants of the phylip format, and support for morphological data sets (i.e., with arbitrary chaaracter state labels). for programs that can work on streams (i.e., one sequence at a time) the reader has been updated. for programs that require the entire alignment to be read in before a task can be performed, the function `std::vector<Sequence> ingest_alignment (std::istream* pios, std::string& alphaName);` (in `seq_reader.h`) can be used. here, `alphaName` is the phyx alphabet type (DNA, AA, BINARY, MULTI), passed by reference, as this is often important information for downstream procedures.
- removal of unnecessary `#include`s throughout
- add missing tests
- phyx now has manpages ([feature request](https://github.com/FePhyFoFum/phyx/issues/56))
    - when writing a new program, use the python script `generate_manpages.py` to generate its manpage. this uses `help2man`, which is easy to install. as long as your `px* -h` output is appropriate (see existing ones) then decent manpages will be geenrated
    - the phyx `Makefile` will (un)install manpages with programs
    - no one will ever use this feature, but it is there!
- all warning/error messages are written to `std::cerr` instead of `std::cout`
- "orphan files" (i.e., those not currently used in the `Makefile` in any way) are moved to the directory `/src/orphan_files`. some of these are deprecated, others are only half implemented (and perhaps forgotten). please resurrect as applicable
- all programs now have a `-C --citation` option. currently these all report the phyx _Bioinformatics_ paper. additional citations should be added for indiviual methods as applicable.
- given the extensive changes, and in preparation for a new release, i've updated the version of the individual programs to `v1.1`

### Issues addressed
- better phylip reader ([bug](https://github.com/FePhyFoFum/phyx/issues/118))
- codon support for `pxclsq` ([feature request](https://github.com/FePhyFoFum/phyx/issues/90))
- fix `pxtrt` leaving a root edge ([bug](https://github.com/FePhyFoFum/phyx/issues/121))
- add manpages ([feature request](https://github.com/FePhyFoFum/phyx/issues/56))
- fix to `pxrmt` speed ([enhancement](https://github.com/FePhyFoFum/phyx/issues/74))
    - this is an issue that appeared in a [paper](https://academic.oup.com/bioinformatics/article/34/6/1053/4582279) for the R package `castor`
    - the issue involves the extensive memory management overhead involved with deleting nodes. `pxrmt` now will use the `pxtrt` tracing procedure is the number of nodes to remove are too large
    - the `castor` example now takes 5 seconds (instead of more than 3 hours; the author actually killed the job at this point)
- `pxrls` now uses regex ([feature request](https://github.com/FePhyFoFum/phyx/issues/59))
- so does `pxrlt` ([feature request](https://github.com/FePhyFoFum/phyx/issues/106))

### New programs
- pxmono: test the monophyletic status of a list of names in a tree (or distribution)
- pxcomp: sequence alignment compositional heterogeneity test (chi-square)
- pxtgen: generate exhaustive (rooted or unrooted) tree topologies for n taxa (limited to 10)
    - for hypothesis testing

### Miscellaneous (no particular order)
- RAII timer added with `timer.h`
    - to use this, do something like:
```c++
include "timer.h"
 ...

void foo (some args) {
    // some initializing code here
    {
        Timer tm;
        // do something here to be timed
    }
    // do something else
}
```
- `pxupgma` was completely rewritten. somewhere along the line, it stopped producing the correct result (using PAUP* as a check)
- `pslssq` now handles all the supported datatypes (DNA, AA, BINARY, MULTI) rather than just the first 2
- `pxboot` overhauled
- `pxclsq` refactored
- `utils.h` now has the following:
```
std::string peek_line (std::istream& pios);
std::vector<std::string> peek_lines (std::istream& pios, const int& n);
```
    - this allows one to peek at line(s) in a stream to make decisions for downstream procedures without losing position in the stream
- `pxaa2cdn` was overhauled (more efficient, and better error-checking)
- more consistent member variable names and formatting
- remedy various `-Wall` compiler warnings
- tree cloner in `tree.h`. important if we want to do something with a tree _and_ keep the original
- `pxs2nex` write the correct Nexus datatype and alphabet
- a bunch of new example files in `/src/TEST`. all the crazy phylip formats, interleaved Nexus files, binary and morphology alignments, etc.
- `pxtrt` overhaul
- about 40K other minor changes
