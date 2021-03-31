# Current programs
Here is a list of the current programs in phyx. If you have suggestions for new programs, please let us know by submitting an [issue](https://github.com/FePhyFoFum/phyx/issues).

By default, all programs are compiled with `make`. An individual program can be compiled with the command `make <name>`.

As with typical Unix programs, help is displayed with the `-h` flag, printing out a menu with the options and their types. Not all of the program options are described below, so make sure to check the help to see what is available. Program version information can be obtained by using the `-V` flag (note uppercase; lowercase `-v` is typically reserved for more verbose output).

A few notes on default behaviour. First, as with standard Unix programs, output files are overwritten without warning, so it is important to be aware of this. Second, most programs that produce output do so to a default format: newick for trees, and fasta for sequences. However, these can always be piped to a subsequent program to change the format (see below).

If you would like to process a large number of files at once a simple command line for loop will let you do that. An example would be if you have a few thousand fasta files that you would like to remove all ambiguous data from, you can run a line of code such as this:
```
for x in *.fa; do pxclsq -s $x -o $x-cln -p 1.0; done
```
The programs can also automatically pipe the output of one into another the input of another. An example of this could be to take the aligned amino acid sequences and guide that to align nucleotide sequences, then clean the file, then create a quick neighbor-joining tree:
```
pxaa2cdn -a amino_acid_alignment -n nucleotide_alignment || pxclsq -p 1.0 || pxnj -o output_tree_file
```

* **pxaa2cdn**: converts AA alignment and unaligned nucleotide to codon alignment

This is a program that lets you change an unaligned nucleotide sequence into its corresponding codon alignment. This is useful for calculating Ka/Ks values and any other analyses that require the nucleotides be aligned by codon. The program requires an aligned amino acid file and the corresponding unaligned nucleotide file. Any sequences that are found in one file and not the other are removed. In addition, the last codon can be removed with the `-r` option.
```
pxaa2cdn -a AA_Alignment.fa -n Unaligned_Nucleotide.fa -o CDN_aln.fa
```

* **pxbdfit**: diversification model inference

Fit a diversification model to an ultrametric tree. Model is controlled by the `-m` flag, which has the options `bd` (default), `yule`, or `best` (the optimal model as determined by AIC). Returns model parameters (b, d, r, e), likelihood, aic, and tree statistics.
```
pxbdfit -t bd.tre -m yule
```

* **pxbdsim**: a birth death simulator

This program is a birth death simulator that allows the user to either specify the number of extant species (`-e`) or the simulation can be run for a given amount of time (`-t`). The simulator also allows for the incorporation of the taxa that have gone extinct (via the `-s` flag).
```
pxbdsim -e 100 -s -b 1 -d 0.5 -o output_tree_file
```

* **pxboot**: sequence alignment resampling (bootstrap or jackknife)

This program is designed to either run a bootstrap or a jackknife on the data. If a bootstrap is run then there is no need to specify an amount but for a jackknife the jackknife percentage must be specified with (`-f`) 

Example jackknife:
```
pxboot -s Alignment -x 112233 -f 0.50 -o output_of_50_jackknife
```
Example bootstrap:
```
pxboot -s Alignment -p parts -o output_of_bootstrap
```
where 'parts' above is a RAxML-style partition file so that conserved-partition bootstrapping can be accomplished.

* **pxbp**: prints out bipartitions that make up the tree

This program takes in a tree file and prints out all the bipartitions that compose the tree.
```
pxbp -t tree_file -o bp_output
```

* **pxbpsq**: prints out bipartions from a sequence file

* **pxcat**: an alignment concatenator

This program is designed to concatenate sequence files and produce a partition file in RAxML/IQtree format. The sequences can be in any format (fasta, phylip, nexus). You can list all the sequences or use a wild card ex ```*.fasta``` which would then concatenate everything that ends in fasta. At times the list of files may be too long and reach the shell limit, in which case you can create a file with all the sequences you could like to be concatenated. An example of this would be if you have a large number of files that end in .fa and you reach a limit from concatenating using *.fa, then you can do ```for x in *.fa;do echo $x >> file.txt;done```. 

```
-s list of sequences, space delimited. Alternatively, a wildcard such as *.fa
-p The name of the parts file to be output
-f a file of all the sequences you'd like concatenated
-u export all characters to uppercase
-o output file (alternatively will output to screen)
```

Concatenating from shell example: ```pxcat -s *.fa *.phy *.nexus -p Parts.model -o Supermatrix.fa```
Concatenating using a file: ```pxcat -f file.txt -p Parts.model -o Supermatrix.fa```

* **pxclsq**: cleanseqs (clean sites based on missing or ambiguous data)

This is a sequence alignment cleaning program, designed to clean based on column occupancy. The proportion you specify with ```-p``` will be the minimum proportion required to be in the column. So if you specify 0.1, then the column will have to have at least 10% occupancy. Assuming a sequence has not more characters, phyx will remove that sequence from the cleaned alignment.

```
-s The name of your sequence file
-o The name of the output file
-p The proportion column occupancy required
-a Specify if your data is amino acids
```

Example where an amino acid alignment needs a minimum of 10% occupancy: ```pxclsq -s Alignment -p 0.1 -a```

* **pxcolt**: collapse poorly-supported edges

This program collapses edges below some support threshold `-l`. This threshold given as a proportion, even if support values are not; e.g., if all edges below a support value of 70 are to be collapsed, specify `-l 0.7`.
```
pxcolt -t tree_file -l threshold
```

* **pxconsq**: a consensus sequence constructor for an alignment
```
pxconsq -s Alignment
```
* **pxcontrates**: a brownian and ou estimator
```

pxcontrates -c contrates_file.txt -t contrates_tree.tre -a 1
```

* **pxfqfilt**: a fastq filter given a mean quality
```

pxfqfilt -s fqfilt_test.fastq -m 10
```

* **pxlog**: a MCMC log manipulator/concatenator

Resamples parameter or tree MCMC samples using some burnin and thinning across an arbitrary number of log files. NOTE: resampling parameters are in terms of number of samples, _not_ number of generations. To determine the attributes of the log files, you can first use the `-i` (`--info`) flag:
```
pxlog -t tree_files -i
```
and then sample accordingly:
```
pxlog -t tree_files -b some_burnin -n some_thinning
```

* **pxlssq**: information about seqs in a file (like ls but for a seq file)

Prints summary information about a sequence alignment, including number of taxa/characters, character frequencies, etc.
```
pxlssq -s Alignment
```
Get statistics for each taxon (rather than the alignment as a whole) by supplying the `-i` argument.

* **pxlstr**: information about trees in a file (like ls but for a tree file)

Prints summary information about trees, including whether it is rooted/binary/ultrametric, number of terminals, tree length, etc.
```
pxlstr -t tree_file
```

* **pxmono**: monophyly tester

Check whether a list of taxa form a monophyletic clade in a tree (or set of trees). The terminal taxa can be specified either with a comma separated list (`-n`) or from a file (one name per line; `-f`).
```
pxmono -t tree_file -n taxon_1,taxon_2,taxon_3
```

* **pxmrca**: information about an mrca

This program will provide the information regarding the most recent common ancestor, giving number of tips in the tree and number of tips for each clade specified. The clade that will be analyzed is the smallest clade containing the tips specified.
```
pxmrca -t mrca_test.tre -m mrca.txt
```

* **pxmrcacut**: a mrca cutter
```
pxmrcacut -t tree -m mrca_file
```

* **pxmrcaname**: a mrca label maker
```
pxmrcaname -t tree -m mrca_file
```

* **pxnj**: Basic neighbor joining program

This is a basic neighbor joining program that will produce trees with the cannonical branch lengths, # of substitutions instead of substitutions per base pair. It allows for parallel processing and is designed to be extremely rapid to give a rough idea of what the final tree will come out to be before doing an ML or Bayesian analysis. For teaching purposes; definitely not for making publishable trees.
```
pxnj -s Alignment.aln
```

* **pxnni**: a nni changer (being trouble shot)

This is a basic nearest neighbor interchange program. Takes in a newick or nexus file with one or more trees and performs a nearest neighbor interchange.
```
pxnni -t tree_file -o output_tree_file
```

* **pxnw**: simple needleman-wunsch

This is a simple alignment program that performs an analysis of pairwise alignments using the Needleman Wunch algorithm for all sequences in the file. It has the options of outputting the scores (`-o`) and the alignment (`-a`). The program also allows the user to input a matrix or uses EDNA for DNA and Blossum62 for AA as defaults.
```
pxnw -s Alignment.aln
```

* **pxpoly**: a polytomy sampler that generates a binary tree.

This program randomly samples 2 descendant lineages from each polytomy in a tree to generate a fully binary tree. Currently restricted to rooted trees.
```
pxpoly -t tree_file
```

* **pxrecode**: a sequence alignment recoder. Currently only to RY-coding, but more coming.

This program will recode a nucleotide alignment using any combination of the recognized recoding schemes: R (A|G), Y (C|T), S (C|G), W (A|T), M (A|C), K (G|T), B (C|G|T), D (A|G|T), H (A|C|T), V (A|C|G). Recoding schemes (e.g., `RY', `SW', `MK', etc.) are specified with the \texttt{-r} argument. If no scheme is provided, RY-coding is used by default.
```
pxrecode -r ry -s Nucleotide.fa
```

* **pxrevcomp**: a reverse complementor
```
pxrevcomp -s Nucleotide.fa
```

* **pxrls**: Taxon relabelling for sequences

Takes two ordered lists of taxon labels, `-c` (current) and `-n` (new), with one label per line. Substitutes the former for the latter in an alignment passed in by `-s` (or stdin). This is convenient for switching between lab codes for analysis and taxon names for data archiving.
```
pxrls -s SeqFile -c CurrentNames -n NewNames
```

* **pxrlt**: Taxon relabelling for trees

Takes two ordered lists of taxon labels, `-c` (current) and `-n` (new), with one label per line. Substitutes the former for the latter in trees passed in by `-t` (or stdin). This is convenient for switching between lab codes for analysis and taxon names for figure preparation.
```
pxrlt -t kingdoms.tre -c kingdoms.oldnames.txt -n kingdoms.newnames.txt
```

* **pxrms**: pruning seqs (like rm but for seqs)

This program is designed to delete sequences from a file, through the input of a file with the sequences you wish to delete that are all on a separate line. Use the `-c` flag to remove the complementary taxa.
```
pxrms -s Nucleotide.fa -f List.txt
```

* **pxrmt**: pruning trees (like rm but for trees)
 
This program is designed to take in a tree and prune tips from the tree that are not wanted, tips can be given either as a comma separated list or can be given in a file. Use the `-c` flag to remove the complementary taxa.
```
pxrmt -t rmt_test.tre -n s1
```

* **pxrr**: rerooting and unrooting trees

This program will re-root trees based off of given outgroup(s) (`-g`) or the program can unroot a tree (`-u`). Outgroups are specified  If not all the outgroups are found in the tree the program will print an error. However, the program can re-root the tree based on the outgroups that _are_ available by using the silent option (`-s`). Alternatively, if the outgroups are ranked in preference but not all necessarily present in a given tree the program can root on the first outgroup present by using the `-r` option. It provides a useful tool for re-rooting thousands of trees which can then be used for analyzing gene discordance across phylogenies.
```
pxrr -t rr_test.tre -g s1,s2
```

* **pxs2fa**: a seq file converter to fasta (force to uppercase with `-u` argument)
```
pxs2fa -s Alignment
```

* **pxs2phy**: a seq file converter to phylip (force to uppercase with `-u` argument)
```
pxs2phy -s Alignment
```

* **pxs2nex**: a seq file converter to nexus (force to uppercase with `-u` argument)
```
pxs2nex -s Alignment
```

* **pxseqgen**: Sequence simulation program

This is a sequence simulator that allows the user to give a tree and specify a model of evolution and sequences will be generated for that tree under the model. Some features are that it allows for the model of evolution to change at nodes along the tree using the (`-m`) option. The program also allows the user to specify rate variation through a value for the shape of the gamma distribution with the (`-g`) option and the user is able to specify the proportion of invariable sites the would like to include using the (`-i`) option. Other options can be found from the help menu by typing (`-h`) after the program.

For multimodel simulations it is easiest to print out the node labels on your tree originally using the (`-p`) option.
Once you know the nodes that you would like the model to change at you can specify these nodes on the input using the (`-m`) option. An example if you wanted two models of evolution on your tree one for the tree and one where it changes at node two, you would enter the command as follows.

if the model you want for the tree is: (.33,.33,.33,.33,.33) where values correspond to (A<->C,A<->G,A<->T,C<->G,C<->T,G<->T)

and the model you want to change to at node two is: (.30,.30,.20,.50,.40) where values correspond to (A<->C,A<->G,A<->T,C<->G,C<->T,G<->T)

The command would be as follows
```
pxseqgen -t tree_file -o output_alignment -m A<->C,A<->G,A<->T,C<->G,C<->T,G<->T,Node#,A<->C,A<->G,A<->T,C<->G,C<->T,G<->T

pxseqgen -t tree_file -o output_alignment -m .33,.33,.33,.33,.33,.33,2,.3,.3,.2,.5,.4,.2
```

* **pxsstat**: multinomial alignment test statistics

This program calculates multinomial alignment test statistics that can be used for assessing model adequacy. Currently limited to the test statistic of Bollback (2002) MBE, but more are coming.
```
pxsstat -s Sequence.fa
```

* **pxstrec**: a state reconstructor

This is a program that does some ancestral state reconstruction and stochastic mapping of categorical characters. There are a number of options and the requirement for a control file. The control file can be as simple as `ancstates = _all_` which designates that you want ancestral states calculated for each node. The can then be output on a tree in a file given by an `-o FILE` option. If you only want to look at particular nodes, these can be designated in the control with the `mrca = MRCANAME tipid1 tipid2`. Then the MRCANAME can be given at the `ancstates = MRCANAME`. If you would like stochastic mapping with the time in the state mapped you can use the same format but instead of `ancstates` you would put `stochtime`. For stochastic number of events `stochnumber` or ``stochnumber_any`. For the stochastic mapping, you will need to designate an MRCA or MRCAs (not _all_). Multiple can be separated by commas or spaces. You can output these to a file with `-n` for number of events, `-a` for the total number of events, and `-m` for the duration. 
```
pxstrec -d test.data.narrow -t test.tre -c config_stmap
```

* **pxsw**: simple smith waterman

This is a simple alignment program that performs an analysis of pairwise alignments using the Smith Waterman algorithm for all sequences in the file. It has the options of outputting the scores (`-o`) and the alignment (`-a`). The program also allows the user to input a matrix or uses EDNA for DNA and Blossum62 for AA as defaults.
```
pxsw -s Alignment.fa
```

* **pxt2new**: a tree file converter to newick
```
pxt2new -t Tree.nex
```

* **pxt2nex**: a tree file converter to vanilla Nexus
```
pxt2nex -t Tree.new
```

* **pxtgen**: exhaustive tree topology generator

This program generates all possible (unrooted or rooted) topologies for `-n` taxa. Limited to <= 10 taxa.
```
pxtgen -n 8
```

* **pxtlate**: Translate nucleotide sequences into amino acids

This program translates nucleotide sequences to their corresponding amino acid sequences. By default it uses the standard translation table, but this can be changed with the `-t` argument (use `-h` to which tables are currently available).
```
pxtlate -s Sequence.fa
```

* **pxtrt**: extract an induced subtree from a larger tree

By specifying a list of taxa, extract the induced subtree from a more comprehensive reference tree. The terminal taxa can be specified either with a comma separated list (`-n`) or from a file (one name per line; `-f`).
```
pxtrt -t tree_file -n taxon_1,taxon_2,taxon_3,taxon_5,taxon_8,taxon_13,taxon_21,taxon_34,taxon_55,taxon_89
```

* **pxtscale**: Tree rescaling.
 
Tree rescaling by providing either scaling factor (`-s`) or root height (`-r`) (not both); the latter requires an
ultrametric tree.
```
pxtscale -t ultra.tre -s 2.0
```

* **pxupgma**: builds a basic upgma tree

This is a basic upgma tree builder, definitely not for making publishable trees, it's here because there are not many out there and it's useful if you're teaching a phylogenetics class and need an example of one of the earliest ways trees were made. It also prints the distance matrix to the screen.
```
pxupgma -s drosophila.aln
```

* **pxvcf2fa**: convert vcf file to fasta alignment (force to uppercase with `-u` argument)
```
pxvcf2fa -s vcf_file
```

# Programs in development
* **pxnprs**: calculate a time calibrated tree with NPRS
* **pxtrsq**: tree-seq; remove from tree and alignment taxa not found in both (JWB)
* **pxgbdb**: a basic genbank database creator
* **pxcomp**: a composition homogeneity test (JWB)
* **pxtdist**: tree distance calculator (JWB)

# Programs in planning
* **pxcoal**: gene trees in species trees; simulation and probabilities. (JWB)
* **pxau**: calculate the "approximately unbiased" test tree probabilities (Shimodaira, 2002). (JWB)