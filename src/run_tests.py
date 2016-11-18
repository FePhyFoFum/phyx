import sys
import os
import subprocess

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def test_program(name):
    x = None
    t = None
    print "TESTING",name
    if name == "pxlstr":
        cm = "pxlstr -t TEST/test.tre"
        t = 'tree #: 0\nrooted: true\nbinary: true\nnterminal: 5\nninternal: 4\nbranch lengths: true\nrttipvar: 0.004634\ntreelength: 1.595\nultrametric: false\nrootheight: NA\n'
    elif name == "pxlssq":
        cm = "pxlssq -s TEST/test.fa"
        t = 'File type: fasta\nNumber of sequences: 5\nIs aligned: true\nSequence length: 20\n--------Nucl TABLE---------\nNucl   Total      Proportion\nA      20         0.2\nC      37         0.37\nG      18         0.18\nT      25         0.25\n-      0          0\nN      0          0\nG+C    55         0.55\n--------Nucl TABLE---------\n'
    elif name == "pxbp":
        cm = "pxbp -t TEST/test.tre"
        t = '1 trees \n3 unique clades found\nTaxonA TaxonB \t1\nTaxonA TaxonB TaxonC \t1\nTaxonD TaxonE \t1\nTSCA: 0\n'
    elif name == "pxbdfit":
        cm = "pxbdfit -t TEST/ultra_100.tre"
        t = 'ntips: 10\nnspeciation: 8\ntreelength: 8.67581\nrootheight: 1.74148\nmodel: bd\nlikelihood: 4.18076\naic: -4.36152\naicc: -2.64723\nb: 1.05759\nd: 0.319273\nr (b-d): 0.738315\ne (d/b): 0.301888\n'
    elif name == "pxtscale":
        cm = "pxtscale -t TEST/test.tre -s 10"
        t = '(((TaxonA:1,TaxonB:0.2999999999999999889):1.25,TaxonC:2.5):4,(TaxonD:2.3000000000000002665,TaxonE:1.6000000000000000888):3);\n'
    elif name == "pxupgma":
        cm = "pxupgma -s TEST/test.fa"
        t = '\tTaxonA\tTaxonB\tTaxonC\tTaxonD\tTaxonE\t\nTaxonA\t0\t18\t12\t14\t14\t\nTaxonB\t18\t0\t16\t16\t14\t\nTaxonC\t12\t16\t0\t15\t15\t\nTaxonD\t14\t16\t15\t0\t9\t\nTaxonE\t14\t14\t15\t9\t0\t\n(((TaxonA:6.000000,TaxonC:6.000000):1.250000,(TaxonD:4.500000,TaxonE:4.500000):2.750000):0.750000,TaxonB:8.000000);\n'
    elif name == "pxs2phy":
        cm = "pxs2phy -s TEST/test.fa"
        t = '5 20\nTaxonA\tAAATTTCCCTGTCCCTTTAA\nTaxonB\tGCTCGAGGGGCCCCAAGACC\nTaxonC\tACGCTCCCCCTTAAAAATGA\nTaxonD\tTCCTTGTTCAACTCCGGTGG\nTaxonE\tTTACTATTCCCCCCCGCCGG\n'
    elif name == "pxfqfilt":
        cm = "pxfqfilt -m 20 -s TEST/test.fastq | wc -l"
        t = "3872\n"
    elif name == "pxs2nex":
        cm = "pxs2nex -s TEST/test.fa"
        t = '#NEXUS\nBEGIN DATA;\n\tDIMENSIONS NTAX=5 NCHAR=20;\n\tFORMAT DATATYPE=DNA INTERLEAVE=NO GAP=-;\n\tMATRIX\n\nTaxonA\tAAATTTCCCTGTCCCTTTAA\nTaxonB\tGCTCGAGGGGCCCCAAGACC\nTaxonC\tACGCTCCCCCTTAAAAATGA\nTaxonD\tTCCTTGTTCAACTCCGGTGG\nTaxonE\tTTACTATTCCCCCCCGCCGG\n;\nend;\n\n'
    elif name == "pxnw":
        cm = "pxnw -s TEST/test.fa | grep TaxonA | grep TaxonB"
        t = 'TaxonA\tTaxonB\t40\n'
    elif name == "pxrecode":
        cm = "pxrecode -s TEST/test.fa"
        t = '>TaxonA\nRRRYYYYYYYRYYYYYYYRR\n>TaxonB\nRYYYRRRRRRYYYYRRRRYY\n>TaxonC\nRYRYYYYYYYYYRRRRRYRR\n>TaxonD\nYYYYYRYYYRRYYYYRRYRR\n>TaxonE\nYYRYYRYYYYYYYYYRYYRR\n'
    elif name == "pxsw":
        cm = "pxsw -s TEST/test.fa | grep TaxonA | grep TaxonB"
        t = 'TaxonA\tTaxonB\t40\n'
    elif name == "pxrms":
        cm = "pxrms -s TEST/test.fa -n TaxonA"
        t = '>TaxonB\nGCTCGAGGGGCCCCAAGACC\n>TaxonC\nACGCTCCCCCTTAAAAATGA\n>TaxonD\nTCCTTGTTCAACTCCGGTGG\n>TaxonE\nTTACTATTCCCCCCCGCCGG\n'
    elif name == "pxboot":
        cm = "pxboot -s TEST/test.fa -x 1"
        t = '>TaxonA\nAAATTCCCCCTGCCCTTTTA\n>TaxonB\nGCTCCGGGGGGCCCCAAGAC\n>TaxonC\nACGCCCCCCCCTAAAAAATA\n>TaxonD\nTCCTTTTTTTAATTCGGGTG\n>TaxonE\nTTACCTTTTTCCCCCGGCCG\n'
    elif name == "pxrevcomp":
        cm = "pxrevcomp -s TEST/test.fa"
        t = '>TaxonA\nTTAAAGGGACAGGGAAATTT\n>TaxonB\nGGTCTTGGGGCCCCTCGAGC\n>TaxonC\nTCATTTTTAAGGGGGAGCGT\n>TaxonD\nCCACCGGAGTTGAACAAGGA\n>TaxonE\nCCGGCGGGGGGGAATAGTAA\n'
    elif name == "pxt2new":
        cm = "pxt2new -t TEST/test_nexus.tre"
        t = '(((TaxonA:0.10000000000000000555,TaxonB:0.02999999999999999889):0.125,TaxonC:0.25):0.4000000000000000222,(TaxonD:0.23000000000000000999,TaxonE:0.16000000000000000333):0.2999999999999999889);\n'
    elif name == "pxbdsim":
        cm = "pxbdsim -e 5 -x 1"
        t = '((taxon_1:0.30855693426109653821,taxon_2:0.30855693426109653821):0.63137104627199058804,(taxon_3:0.93851839959563365667,(taxon_4:0.53977584969885905597,taxon_5:0.53977584969885905597):0.39874254989677465622):0.0014095809374534387821);\n'
    elif name == "pxs2fa":
        cm = "pxs2fa -s TEST/Concat_Sequence2.NEX"
        t = '>Sequence1\nAAATTTCCCTTTCCCTTTAAA\n>Sequence2\nGGGGGGGGGGCCCCCCCCCCA\n>Sequence3\nCCCCCCCCCCCCAAAAAAAAA\n>Sequence9\nAAATTTCCCTTTCCCTTTAAA\n>Sequence10\nGGGGGGGGGGCCCCCCCCCCA\n>Sequence11\nCCCCCCCCCCCCAAAAAAAAA\n>Sequence8\nTTTTTTTTCCCCCCCGGGGGA\n'
    elif name == "pxrr":
        cm = "pxrr -t TEST/test.tre -g TaxonA"
        t = '(TaxonA:0.050000000000000002776,(TaxonB:0.02999999999999999889,((TaxonD:0.23000000000000000999,TaxonE:0.16000000000000000333):0.69999999999999995559,TaxonC:0.25):0.125):0.050000000000000002776);\n'
    elif name == "pxnj":
        cm = "pxnj -s TEST/test.fa"
        t = '((((TaxonA:0.300000,TaxonC:0.300000):0.112500,TaxonB:0.437500):0.087500,TaxonD:0.250000):0.100000,TaxonE:0.100000);\n'
    elif name == "pxconsq":
        cm = "pxconsq -s TEST/test.fa"
        t = '>consensus\nRHNYKNBBSNNYHMMRNHVV\n'
    elif name == "pxseqgen":
        cm = "pxseqgen -t TEST/test.tre -x 1 -l 1"
        t = '>TaxonE\nT\n>TaxonD\nT\n>TaxonC\nA\n>TaxonB\nT\n>TaxonA\nA\n'
    elif name == "pxrmt":
        cm = "pxrmt -t TEST/test.tre -n TaxonA"
        t = '((TaxonC:0.25,TaxonB:0.15499999999999999889):0.4000000000000000222,(TaxonD:0.23000000000000000999,TaxonE:0.16000000000000000333):0.2999999999999999889);\n'
    elif name == "pxcat":
        cm = "pxcat -s TEST/test.fa TEST/test.fa"
        t = '>TaxonA\nAAATTTCCCTGTCCCTTTAAAAATTTCCCTGTCCCTTTAA\n>TaxonB\nGCTCGAGGGGCCCCAAGACCGCTCGAGGGGCCCCAAGACC\n>TaxonC\nACGCTCCCCCTTAAAAATGAACGCTCCCCCTTAAAAATGA\n>TaxonD\nTCCTTGTTCAACTCCGGTGGTCCTTGTTCAACTCCGGTGG\n>TaxonE\nTTACTATTCCCCCCCGCCGGTTACTATTCCCCCCCGCCGG\n'
    else:
        return
    p = subprocess.Popen(cm,shell=True,stdout=subprocess.PIPE)
    x = p.communicate()
    #print x
    x = x[0]
    if x == t:
        print bcolors.OKBLUE+"PASSED"+bcolors.ENDC
        return True
    else:
        print bcolors.FAIL+"FAILED"+bcolors.ENDC
        return False

if __name__ == "__main__":
    if len(sys.argv) != 1:
        print "python run_tests.py"
        sys.exit(0)
    
    passed = 0
    failed = 0
    failedl = []
    print "================="
    for i in os.listdir("."):
        if i[:2] == "px":
            t = test_program(i)
            if t == True:
                passed += 1
            elif t == False:
                failed += 1
                failedl.append(i)
            else:
                print bcolors.WARNING+"no test for :"+i+bcolors.ENDC
            print "================="
    print "PASSED TESTS:",passed
    print "FAILED TESTS:",failed
    if failed > 0:
        print "These failed:"
        print "\t",",".join(failedl)
