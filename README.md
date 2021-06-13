[![Build Status](https://travis-ci.com/FePhyFoFum/phyx.svg?branch=master)](https://travis-ci.com/FePhyFoFum/phyx)

<img src="logo.png" alt="phyx logo" width=208px/> 

---

**_Note_** Phyx recently overwent an overhaul such that a simple `git pull && git make` will fail. Instead, see instructions [here](#problems-after-updating-git-pull).

---

**phyx** performs phylogenetics analyses on trees and sequences. See installation instructions for Linux and Mac including any dependencies on the wiki [here](https://github.com/FePhyFoFum/phyx/wiki/Installation) or below.

Authors: Joseph W. Brown\*, Joseph F. Walker\*, and Stephen A. Smith (\* equal contribution)

Citation: [Brown, J. W., J. F. Walker, and S. A. Smith; Phyx: phylogenetic tools for unix. Bioinformatics 2017; 33 (12): 1886-1888. doi: 10.1093/bioinformatics/btx063](https://academic.oup.com/bioinformatics/article/33/12/1886/2975328/Phyx-phylogenetic-tools-for-unix)

License: GPL https://www.gnu.org/licenses/gpl-3.0.html

Some of the sequence comparison operations use the very nice [edlib library](https://github.com/Martinsos/edlib#alignment-methods). These are reported in this publication: [Martin Šošić, Mile Šikić; Edlib: a C/C ++ library for fast, exact sequence alignment using edit distance. Bioinformatics 2017 btw753. doi: 10.1093/bioinformatics/btw753](https://academic.oup.com/bioinformatics/article/33/9/1394/2964763/Edlib-a-C-C-library-for-fast-exact-sequence).

## Documentation
Documentation resides in several locations (all slightly out of date, alas). A [pdf manual](https://github.com/FePhyFoFum/phyx/tree/master/doc) is available in the `doc/` directory. A slightly-less-out-of-date list of the current programs with examples can be found [on the wiki](https://github.com/FePhyFoFum/phyx/wiki/Program-list). Help for individual programs can be obtained with either `PROGRAM -h` or (if installed, see below) `man PROGRAM`. See a brief overview [here](https://twitter.com/i/moments/1067839564927008769).

Still, there are a whack of programs, so it can be difficult to remember the name of which program does what. Here is a quick reference:

Program | Short description
:--- | :---
pxaa2cdn | produce a codon alignment from and AA alignment and unaligned nucleotides
pxbdfit | diversification model inference
pxbdsim | a birth death simulator
pxboot | sequence alignment resampling (bootstrap or jackknife)
pxbp | prints out bipartitions that make up the tree
pxcat | an alignment concatenator
pxclsq | clean sites based on missing or ambiguous data
pxcltr | general tree cleaner
pxcolt | collapse poorly-supported edges
pxcomp | a composition homogeneity test
pxcomp | compositional homogeneity test
pxconsq | a consensus sequence constructor for an alignment
pxcontrates | a brownian and ou estimator
pxfqfilt | a fastq filter given a mean quality
pxlog | a MCMC log manipulator/concatenator
pxlssq | information about seqs in a file (like ls but for an alignment file)
pxlstr | information about trees in a file (like ls but for a tree file)
pxmono | monophyly tester
pxmrca | information about an mrca
pxmrcacut | a mrca cutter
pxmrcaname | a mrca label maker
pxnj | neighbour-joining tree inference
pxnni | a nni changer
pxnw | needleman-wunsch alignment
pxpoly | a polytomy sampler that generates a binary tree
pxrecode | a sequence alignment recoder
pxrevcomp | a reverse complementor
pxrls | taxon relabelling for sequences
pxrlt | taxon relabelling for trees
pxrmk | remove two-degree nodes from a tree
pxrms | pruning seqs (like rm but for seqs)
pxrmt | pruning trees (like rm but for trees)
pxrr | rerooting and unrooting trees
pxs2fa | convert an alignment to fasta format
pxs2nex | convert an alignment to nexus format
pxs2phy | convert an alignment to phylip format
pxseqgen | Sequence simulation program
pxssort | sequence sorter
pxsstat | multinomial alignment test statistics
pxstrec | a state reconstructor
pxsw | smith waterman alignment
pxt2new | convert a tree to newick format
pxt2nex | convert a tree to vanilla Nexus format
pxtcol | annotate tree to colour edges
pxtcomb | tree combiner
pxtgen | exhaustive tree topology generator
pxtlate | translate nucleotide sequences into amino acids
pxtrt | extract an induced subtree from a larger tree
pxtscale | tree rescaling
pxupgma | upgma tree inference
pxvcf2fa | convert vcf file to fasta alignment

## Problems after updating (git pull)
If you have been using phyx and things are not working after a recent pull, this is because of a change in configuration. Please do the following in the `src` directory to remedy the situation:

    make distclean
    autoreconf -fi
    ./configure
    make
    make check
    sudo make install

# Installation instructions 
phyx requires a few dependencies. Since installation of these dependencies differs on [Linux](#linux-install) vs. [Mac OSX](#mac-install), we've separated the instructions below. 

## Mac install
Mac has become increasingly difficult to support at the command line with changes every version on location and standards for compilation tools. First, distribution of compiled programs is very difficult. Furthermore, Mac now defaults to clang as a C/C++ compiler, which does not support OpenMP.  For **Mac OSX 10.12**, we have found that you can install with clang using the simple instructions and [homebrew](http://brew.sh/) *or* using a fresh installation of gcc from [here](http://hpc.sourceforge.net/). Instructions for both are below (don't use both, choose one, probably the simple one). For simple instructions click [here](#binary-install-with-homebrew), and for advanced instructions click [here](#install-with-hpc-gcc-advanced-instructions).

### Binary install with Homebrew

1. Install the Homebrew package manager:

        /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

2. Install the Brewsci phyx package:

        brew install brewsci/bio/phyx

### Build from source with Homebrew

1. Install the Homebrew package manager:

        /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

2. Install dependencies from homebrew:

        brew install git cmake nlopt armadillo

3. On to phyx. first, clone the repository (if you haven't already):

        git clone https://github.com/FePhyFoFum/phyx.git

4. Install phyx

        cd phyx/src
        autoconf
        ./configure
        make
        make check

If you want to install it so it is available anywhere in your system, do:

        sudo make install

### Install with HPC GCC (advanced instructions)
1. Install gcc and gfortran. Download gcc-6.2-bin.tar.gz or more recent from http://hpc.sourceforge.net/. Install with:
    ​     

        sudo tar -xvf gcc-6.2-bin.tar -C /

2. Install autoconf from http://ftp.gnu.org/gnu/autoconf/. Get autoconf-latest.tar.gz, then:

        tar -xzf autoconf-latest.tar.gz
        cd autoconf-2.69
        ./configure --prefix=/usr/local/autoconf-2.69
        make
        sudo make install
        ln -s autoconf-2.69 /usr/local/autoconf
    
3. On to phyx. first, clone the repository (if you haven't already):

        git clone https://github.com/FePhyFoFum/phyx.git

4. Install cmake and install Armadillo. Get cmake from https://cmake.org/download/. I got https://cmake.org/files/v3.6/cmake-3.6.2-Darwin-x86_64.tar.gz. Get armadillo from the `deps` directory or http://arma.sourceforge.net/download.html, get the stable one. Untar it. Double click the Cmake.app. Click "Browse source..." and choose the armadillo folder that was created after untaring. Click "Browse build..." and choose the same folder as browse source. Click "Configure" and then click "Generate". Go to the terminal and browse to that armadillo folder and type:

        make
        sudo make install

5. Install nlopt. Get armadillo from the `deps` directory or go to http://ab-initio.mit.edu/wiki/index.php/NLopt#Download_and_installation and download the latest (probably nlopt-2.4.2.tar.gz). Untar and browse in the terminal to that directory:

        ./configure --without-octave --without-matlab
        make
        sudo make install

6. Compile phyx. Now you can go to the src directory of phyx and type:

        autoconf
        ./configure
        make
        make check
        sudo make install

and all the programs should compile without issue. 


## Linux install

These instructions work for most ubuntu versions as well as debian. 

1. Install general dependencies:

        sudo apt-get install git autotools-dev autoconf automake cmake libtool liblapack-dev libatlas-cpp-0.6-dev libnlopt-cxx-dev

2. Clone the phyx repo (if you haven't already):

        git clone https://github.com/FePhyFoFum/phyx.git

3. Install armadillo dependency  

**Note**: it is possible to get from apt-get, but need version >= 5.2:

        sudo apt-get install libarmadillo-dev

On debian it was necessary to use backports:

        sudo apt-get -t jessie-backports install libarmadillo-dev

If that is not possible, compile the provided code:

        cd phyx/deps
        tar -xvzf armadillo-7.400.2.tgz
        cd armadillo-7.400.2
        ./configure
        make
        sudo make install

4. Finally, install phyx:

        cd phyx/src
        autoconf
        ./configure
        make
        make check

If you want to install it so it is available anywhere in your system, do:

        sudo make install

