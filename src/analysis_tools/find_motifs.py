import os, sys

from collections import namedtuple, OrderedDict, defaultdict

import re

import numpy

import pysam

try: import grit
except ImportError: sys.path.insert(0, "/home/nboley/grit/grit/")

from grit.files.reads import CAGEReads, RAMPAGEReads

VERBOSE = False
QUIET = False
FIX_CHRM_NAMES_FOR_UCSC = False
NTHREADS = 1

TSSLoc = namedtuple('TSS', ['chrm', 'strand', 'start', 'stop'])

FLANK_SIZE = 50

motif_names = """
BREu
TATA1
TATA2
BREd
XCPE1
MTE
Inr
DCE
DPE""".strip().split()

motif_pats = """
[GC]GC[GA]GGCC
[CG]TATA[AT]A[AT][AG]
T[AG]G[CT][ACGT][ACGT]AGTGG
[CG][AG][AG]CGCC
[GAT][GC]G[TC]GG[GA]A[GC][AC]
C[CG]A[AG]C[CG][CG]AAC
[CT][CT]A[ACGT][AT][CT][CT]
CTTC.{3,40}?CTGT.{3,40}?AGC
[GT]CGGTT[CG][GT]""".strip().split()

motifs = OrderedDict(zip(motif_names, (re.compile(pat) for pat in motif_pats)))

def search_for_motifs(seq):
    matches = defaultdict(list)
    for motif, pattern in motifs.items():
        for match in re.finditer(pattern, seq):
            matches[motif].append(match.span()[0])
    return matches

try: rev_comp_table = str.maketrans("ACGT", "TGCA")
except:
    import string
    rev_comp_table = string.maketrans("ACGT", "TGCA")

class TSS(TSSLoc):
    def get_flanking_seq(self, fasta):
        pos = self.start + (self.stop-self.start)/2
        region = (self.chrm, pos-FLANK_SIZE, pos+FLANK_SIZE)
        seq = str(fasta.fetch(*region).upper())
        if self.strand == '-':
            return seq.translate(rev_comp_table)[::-1] 
        else:
            return seq

def load_TSSs(fp):
    tss_s = []
    for line in fp:
        if line.startswith("track") or line.startswith("#"): continue
        data = line.split()
        assert data[5] in '+-', "Invlaid strand '{}'".format(data[5])
        tss_s.append( TSS(data[0], data[5], int(data[1]), int(data[2])) )
    return tss_s

def log(msg, level=None):
    if QUIET: return
    if level == None or (level == 'VERBOSE' and VERBOSE):
        print >> sys.stderr, msg

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Annotate promoters associated with experimental TSSs.')

    parser.add_argument( 'TSSs', type=argparse.FileType("r"),
        help='Narrow peak file containing annotated TSSs.')

    parser.add_argument( 'fasta', type=pysam.Fastafile,
        help='Indexed fasta file with the reference sequence.')

    parser.add_argument( '--cage-reads', type=argparse.FileType('rb'), 
        help='BAM file containing mapped cage reads.')
    parser.add_argument( '--cage-read-type', 
                         choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the reads that maps to the genome without being reverse complemented are assumed to be on the '+'. default: auto")

    parser.add_argument( '--rampage-reads', type=argparse.FileType('rb'), 
        help='BAM file containing mapped rampage reads.')
    parser.add_argument( '--rampage-read-type', 
                         choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the first read in a pair that maps to the genome without being reverse complemented are assumed to be on the '+' strand. default: auto")
    
    parser.add_argument( '--out-fname', '-o', 
                         help='Output file name. (default stdout)')
    
    parser.add_argument( '--ucsc', default=False, action='store_true', 
        help='Format the contig names to work with the UCSC genome browser.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', 
                         help='Whether or not to print status information.')
    parser.add_argument( '--quiet', '-q', default=False, action='store_true', 
                         help='Suppress all messages.')

    parser.add_argument( '--threads', '-t', default=1, type=int,
                         help='The number of threads to run.')

        
    args = parser.parse_args()

    global VERBOSE
    if args.verbose: VERBOSE = True 

    global QUIET
    if args.quiet: 
        QUIET = True 
        VERBOSE = False
    
    global FIX_CHRM_NAMES_FOR_UCSC
    FIX_CHRM_NAMES_FOR_UCSC = args.ucsc
    
    global NTHREADS
    NTHREADS = args.threads
        
    if args.cage_reads != None:
        assert args.rampage_reads == None, \
            "Can not use both RAMPAGE and CAGE reads"
        if VERBOSE: 
            print( "Loading %s" % args.cage_reads.name )
        rev_reads = {'forward':False, 'backward':True, 'auto': None}[
            args.cage_read_type]
        promoter_reads = CAGEReads(args.cage_reads.name, "rb").init(
            reverse_read_strand=rev_reads, ref_genes=ref_genes)
    elif args.rampage_reads != None:
        assert args.cage_reads == None, "Can not use RAMPAGE and CAGE reads"
        if VERBOSE: 
            log( "Loading %s" % args.rampage_reads.name )
        rev_reads = {'forward':False, 'backward':True, 'auto': None}[
            args.rampage_read_type]
        promoter_reads = RAMPAGEReads(args.rampage_reads.name, "rb").init(
            reverse_read_strand=rev_reads)
    else:
        promoter_reads = None
        pass
    
    output_stream = ( open(args.out_fname, "w") 
                      if args.out_fname != None
                      else sys.stdout )
    
    return ( args.TSSs, args.fasta,
             promoter_reads, output_stream )

def main():
    tss_s_fp, fasta, reads, ofp = parse_arguments()
    tss_s = load_TSSs(tss_s_fp)
    """
    for tss in tss_s:
        cov = reads.build_read_coverage_array(
            tss.chrm, tss.strand, tss.start, tss.stop)
        print tss, cov
        return
    """
    rv = OrderedDict( (motif, numpy.zeros([2*FLANK_SIZE], dtype=float)) 
                      for motif in motifs.keys() )

    for i, tss in enumerate(tss_s):
        if i > 0 and i%1000 == 0: 
            print("Finished {}/{}".format(i, len(tss_s)) )
        seq = tss.get_flanking_seq(fasta)
        matches = search_for_motifs( seq )
        for motif, hits in matches.items():
            for hit in hits:
                rv[motif][hit] += 1./len(hits)

    for motif, cnts in rv.items():
        with open("Frac.{}.txt".format(motif), "w") as ofp:
            for cnt in cnts:
                print >> ofp, "{:e}".format(float(cnt)/len(cnts))
    # load the peaks list
    
    # load the fasta file
    
    # load the bam(s)
    pass

if __name__ == '__main__':
    main()
