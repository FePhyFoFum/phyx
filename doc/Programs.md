# Current programs
Here is a list of the current programs in phyx. If you have suggestions for new programs, please let us know by submitting an [issue](https://github.com/FePhyFoFum/phyx/issues).

By default, all programs are compiled with `make`. An individual program can be compiled with the command `make <nprogram ame>`.

As with typical Unix programs, help is displayed with the `-h` flag, printing out a menu with the options and their types. Alternatively, help information can be displayed in individual man pages `man <program name>`. Not all of the program options are described below, so make sure to check the help to see what is available. Program version information can be obtained by using the `-V` flag (note uppercase; lowercase `-v` is typically reserved for more verbose output).

A few notes on default behaviour. First, as with standard Unix programs, output files are overwritten without warning, so it is important to be aware of this. Second, most programs that produce output do so to a default format: newick for trees, and fasta for sequences. However, these can always be piped to a subsequent program to change the format (see below). By default, all results are written to standard output (i.e., the terminal screen), as this is what allows the programs to be piped together. However, all programs have the option of writing to a file with the `-o` option. 

If you would like to process a large number of files at once a simple command line for loop will let you do that. An example would be if you have a few thousand fasta files that you would like to remove all ambiguous data from, you can run a line of code such as this:
```
for x in *.fa; do pxclsq -s $x -o $x-cln -p 1.0; done
```
The programs can also automatically pipe the output of one into another the input of another. An example of this could be to take the aligned amino acid sequences and guide that to align nucleotide sequences, then clean the file, then create a quick neighbour-joining tree:
```
pxaa2cdn -a amino_acid_alignment -n nucleotide_alignment || pxclsq -p 1.0 || pxnj -o output_tree_file
```

* **pxaa2cdn**: produce a codon alignment from and AA alignment and unaligned nucleotides

This is a program that lets you change an unaligned nucleotide sequence into its corresponding codon alignment. This is useful for calculating Ka/Ks values and any other analyses that require the nucleotides be aligned by codon. The program requires an aligned amino acid file and the corresponding unaligned nucleotide file. Any sequences that are found in one file and not the other are removed. In addition, the last codon can be removed with the `-r` option.
```
pxaa2cdn -a aa_alignment_file -n nucleotide_alignment_file
```

* **pxbdfit**: diversification model inference

Fit a diversification model to an ultrametric tree. Model is controlled by the `-m` flag, which has the options `bd` (default), `yule`, or `best` (the optimal model as determined by AIC). Returns model parameters (b, d, r, e), likelihood, aic, and tree statistics.
```
pxbdfit -t tree_file -m yule
```

* **pxbdsim**: a birth death simulator

This program is a birth death simulator that allows the user to either specify the number of extant species (`-e`) or the simulation can be run for a given amount of time (`-t`). The simulator also allows for the incorporation of the taxa that have gone extinct (via the `-s` flag).
```
pxbdsim -e 100 -s -b 1 -d 0.5
```

* **pxboot**: sequence alignment resampling (bootstrap or jackknife)

This program is designed to either run a bootstrap (sampling with replacement; default) or a jackknife (sampling without replacement) on an alignment. To perform a jackknife a percentage must be specified with the `-f` argument.

Example jackknife:
```
pxboot -s alignment_file -x 112233 -f 0.50
```
Example bootstrap:
```
pxboot -s alignment_file -p partition_file
```
where 'partition_file' above is a RAxML-style partition file so that conserved-partition bootstrapping can be accomplished.

* **pxbp**: prints out bipartitions that make up the tree

This program takes in a tree file and prints out all the bipartitions that compose the tree. See the help for all program options.
```
pxbp -t tree_file
```

* **pxcat**: an alignment concatenator

Concatenates sequence files and produces a partition file in RAxML/IQtree format. The sequences can be in any format (fasta, phylip, nexus). You can list all the sequences or use a wild card ex `*.fasta` which would then concatenate everything that ends in fasta. At times the list of files may be too long and reach the shell limit, in which case you can create a file with all the sequences you could like to be concatenated. An example of this would be if you have a large number of files that end in .fa and you reach a limit from concatenating using *.fa, then you can do:
```
for x in *.fa;do echo $x >> file.txt;done
```


Concatenating from shell:
```
pxcat -s *.fa *.phy *.nexus -p output_partition_file -o concatenated_alignment
```
Concatenating using a file:
```
pxcat -f file.txt -p output_partition_file -o concatenated_alignment
```

* **pxclsq**: clean sites based on missing or ambiguous data

Remove alignment sites based on proportion of missing data (`-p`). By default, removes individual sites from all sequences. Alternatively, missing data can be considered by codon (`-c`) or by taxa (`-t`, where entire sequences with missing data > p are removed).

```
pxclsq -s alignment_file -p 0.1
```

* **pxcltr**: general tree cleaner

Removes annotations (node labels, `-l`), 'knuckles' (2-degree nodes; `-k`), and root edges (`-r`) to generate a 'vanilla' newick representation. By default removes all properties. Alternatively choose 1 property.
```
pxcltr -t tree_file
```

* **pxcolt**: collapse poorly-supported edges

This program collapses edges below some support threshold `-l`. This threshold given as a proportion, even if support values are not; e.g., if all edges below a support value of 70 are to be collapsed, specify `-l 0.7`.
```
pxcolt -t tree_file -l threshold
```

* **pxcomp**: compositional homogeneity test

Chi-square test for equivalent character state counts across lineages.
```
pxcomp -s alignment_file
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
pxfqfilt -s fastq_file -m 10
```

* **pxlog**: a MCMC log manipulator/concatenator

Resamples parameter or tree MCMC samples using some burnin and thinning across an arbitrary number of log files. NOTE: resampling parameters are in terms of number of samples, _not_ number of generations. Can combine burnin (`-b`) and thinning (`-n`). To determine the attributes of the log files, you can first use the `-i` (`--info`) flag:
```
pxlog -t tree_log_file -i
```
and then sample accordingly:
```
pxlog -t tree_log_file -b 10000 -n 100
```

* **pxlssq**: information about seqs in a file (like ls but for an alignment file)

Prints summary information about a sequence alignment, including number of taxa/characters, character frequencies, etc. Get statistics for each taxon (rather than the alignment as a whole) by supplying the `-i` argument.
```
pxlssq -s alignment_file
```

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
pxmrca -t tree_file -m mrca.txt
```

* **pxmrcacut**: a mrca cutter

Extract a subclade from a tree through specification of an MRCA.
```
pxmrcacut -t tree_file -m mrca_file
```

* **pxmrcaname**: a mrca label maker

Label internal nodes (specified by MRCA statements) with clade names from a file (`-m`).
```
pxmrcaname -t tree_file -m mrca_file
```

* **pxnj**: neighbour-joining tree inference

A basic neighbour-joining program that will produce trees with the 'canonical' branch lengths (# of substitutions instead of average # of substitutions per site). For teaching purposes; definitely not for making publishable trees.
```
pxnj -s alignment_file
```

* **pxnni**: a nni changer

Performs a nearest neighboor interchange (NNI) on a tree.
```
pxnni -t tree_file
```

* **pxnw**: needleman-wunsch alignment

An alignment program that performs an analysis of pairwise alignments using the Needleman Wunch algorithm for all sequences in the file. It has the options of outputting the scores (`-o`) and the alignment (`-a`). The program also allows the user to input a matrix or uses EDNA for DNA and Blossum62 for AA as defaults.
```
pxnw -s alignment_file
```

* **pxpoly**: a polytomy sampler that generates a binary tree

Randomly samples 2 descendant lineages from each polytomy in a tree to generate a fully binary tree. Currently restricted to rooted trees.
```
pxpoly -t tree_file
```

* **pxrecode**: a sequence alignment recoder

Recodes a nucleotide alignment using any combination of the recognized recoding schemes: R (A|G), Y (C|T), S (C|G), W (A|T), M (A|C), K (G|T), B (C|G|T), D (A|G|T), H (A|C|T), V (A|C|G). Recoding schemes (e.g., `RY`, `SW`, `MK`, etc.) are specified with the `-r` argument. If no scheme is provided, RY-coding is used by default.
```
pxrecode -r ry -s alignment_file
```

* **pxrevcomp**: a reverse complementor

Replace sequences with their reverse complement. By default this is applied to all sequences; individual sequences can be specified with the `-i` argument.
```
pxrevcomp -s alignment_file
```

* **pxrls**: taxon relabelling for sequences

Takes two ordered lists of taxon labels, `-c` (current) and `-n` (new), with one label per line. Substitutes the former for the latter in an alignment passed in by `-s` (or stdin). This is convenient for switching between lab codes for analysis and taxon names for data archiving.
```
pxrls -s alignment_file -c curren_name_file -n new_name_file
```

* **pxrlt**: taxon relabelling for trees

Takes two ordered lists of taxon labels, `-c` (current) and `-n` (new), with one label per line. Substitutes the former for the latter in trees passed in by `-t` (or stdin). This is convenient for switching between lab codes for analysis and taxon names for figure preparation.
```
pxrlt -t tree_file -c curren_name_file -n new_name_file
```

* **pxrmk**: remove two-degree nodes from a tree

Removes two-degree internal nodes ('knuckles') from a tree.
```
pxrmk -t tree_file
```

* **pxrms**: pruning seqs (like rm but for seqs)

Delete all sequences specified by a comma-deleimted list (`-n`) or a file (`-f`; one per line) from an alignment. Use the `-c` flag to remove the complementary taxa.
```
pxrms -s alignment_file -f taxon_name_file
```

* **pxrmt**: pruning trees (like rm but for trees)

Prune all taxa specified by a comma-deleimted list (`-n`) or a file (`-f`; one per line) from a tree. Use the `-c` flag to remove the complementary taxa.
```
pxrmt -t tree_file -n s1
```

* **pxrr**: rerooting and unrooting trees

Re-root trees based off of given outgroup(s) (`-g`) or unroot (`-u`). Outgroups are specified  If not all the outgroups are found in the tree the program will print an error. However, the program can re-root the tree based on the outgroups that _are_ available by using the silent option (`-s`). Alternatively, if the outgroups are ranked in preference but not all necessarily present in a given tree the program can root on the first outgroup present by using the `-r` option.
```
pxrr -t tree_file -g s1,s2
```

* **pxs2fa**: convert an alignment to fasta format

Convert characters to uppercase with the `-u` argument.
```
pxs2fa -s alignment_file
```

* **pxs2nex**: convert an alignment to nexus format

Convert characters to uppercase with the `-u` argument.
```
pxs2nex -s alignment_file
```

* **pxs2phy**: convert an alignment to phylip format

Convert characters to uppercase with the `-u` argument.
```
pxs2phy -s alignment_file
```

* **pxseqgen**: Sequence simulation program

This is a sequence simulator that allows the user to give a tree and specify a model of evolution and sequences will be generated for that tree under the model. Some features are that it allows for the model of evolution to change at nodes along the tree using the (`-m`) option. The program also allows the user to specify rate variation through a value for the shape of the gamma distribution with the (`-g`) option and the user is able to specify the proportion of invariable sites the would like to include using the (`-i`) option. A rate matrix can be specified with the `-r` option () in the order: A<->C,A<->G,A<->T,C<->G,C<->T,G<->T). Other options can be found from the help menu by typing (`-h`) after the program.

For multimodel simulations it is easiest to print out the node labels on your tree originally using the (`-p`) option. Once you know the nodes that you would like the model to change at you can specify these nodes on the input using the (`-m`) option with the format: rates1,node#,rates2. For example, if the model you want for the tree is (.33,.33,.33,.33,.33,.33) and the model you want to change to at node two is (.3,.3,.2,.5,.4,.2):
```
pxseqgen -t tree_file -m .33,.33,.33,.33,.33,.33,2,.3,.3,.2,.5,.4,.2
```

* **pxssort**: sequence sorter

Sort sequences by id (`-b 1`; default), reverse id (`-b 2`), length (`-b 3`), or reverse length (`-b 4`).
```
pxssort -s alignment_file -b 3
```

* **pxsstat**: multinomial alignment test statistics

This program calculates multinomial alignment test statistics that can be used for assessing model adequacy. Currently limited to the test statistic of Bollback (2002) MBE, but more are coming.
```
pxsstat -s alignment_file
```

* **pxstrec**: a state reconstructor

Ancestral state reconstruction and stochastic mapping of categorical characters. There are a number of options and the requirement for a control file. The control file can be as simple as `ancstates = _all_` which designates that you want ancestral states calculated for each node. The can then be output on a tree in a file given by an `-o FILE` option. If you only want to look at particular nodes, these can be designated in the control with the `mrca = MRCANAME tipid1 tipid2`. Then the MRCANAME can be given at the `ancstates = MRCANAME`. If you would like stochastic mapping with the time in the state mapped you can use the same format but instead of `ancstates` you would put `stochtime`. For stochastic number of events `stochnumber` or ``stochnumber_any`. For the stochastic mapping, you will need to designate an MRCA or MRCAs (not _all_). Multiple can be separated by commas or spaces. You can output these to a file with `-n` for number of events, `-a` for the total number of events, and `-m` for the duration. 
```
pxstrec -d test.data.narrow -t tree_file -c config_stmap
```

* **pxsw**: smith waterman alignment

An alignment program that performs an analysis of pairwise alignments using the Smith Waterman algorithm for all sequences in the file. It has the options of outputting the scores (`-o`) and the alignment (`-a`). The program also allows the user to input a matrix or uses EDNA for DNA and Blossum62 for AA as defaults.
```
pxsw -s alignment_file
```

* **pxt2new**: convert a tree to newick format
```
pxt2new -t tree_file
```

* **pxt2nex**: convert a tree to vanilla Nexus format
```
pxt2nex -t tree_file
```

* **pxtcol**: annotate tree to colour edges

Uses tab-delimited files of mrcas (`-m`) or nodeids (`-d`) to annotate a tree string with colouring information. Output is a Nexus-formatted tree file which can be read by FigTree.
```
pxtol -t tree_file -d node_annotation_file.tsv
```

* **pxtcomb**: tree combiner

Combine a set of trees from one file (`-a`) into a tree from another (`-t`).
```
pxtcomb -t tree_file -a alternate_tree_file
```

* **pxtgen**: exhaustive tree topology generator

Generate all possible unrooted (default) or rooted (`-r`) topologies for `-n` taxa. Limited to <= 10 taxa. Report just the number of possible trees with the `-c` option.
```
pxtgen -n 8 -r
```

* **pxtlate**: Translate nucleotide sequences into amino acids

Translates nucleotide sequences to their corresponding amino acid sequences. By default it uses the standard translation table, but this can be changed with the `-t` argument (use `-h` to which tables are currently available).
```
pxtlate -s alignment_file
```

* **pxtrt**: extract an induced subtree from a larger tree

By specifying a list of taxa, extract the induced subtree from a more comprehensive reference tree. The terminal taxa can be specified either with a comma separated list (`-n`) or from a file (one name per line; `-f`).
```
pxtrt -t tree_file -n taxon_1,taxon_2,taxon_3,taxon_5,taxon_8,taxon_13,taxon_21,taxon_34,taxon_55,taxon_89
```

* **pxtscale**: Tree rescaling.

Tree rescaling by providing either scaling factor (`-s`) or root height (`-r`) (not both); the latter requires an ultrametric tree.
```
pxtscale -t tree_file -s 2.0
```

* **pxupgma**: upgma tree inference

This is a basic upgma tree builder based on uncorrected p-distances. Definitely not for making publishable trees; it's here because there are not many out there and it's useful if you're teaching a phylogenetics class and need an example of one of the earliest ways trees were made.
```
pxupgma -s alignment_file
```

* **pxvcf2fa**: convert vcf file to fasta alignment

Convert characters to uppercase with the `-u` argument.
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