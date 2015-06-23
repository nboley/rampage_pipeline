import os, sys

from collections import namedtuple, OrderedDict, defaultdict

import re

import numpy

import pysam

sys.path.insert(0, "/users/nboley/src/TF_binding/")
#from motif_tools import Motif, iter_motifs
from load_motifs import load_motifs 

motif_names = """
TATA
BDP1
HMGN3
GTF2I
GTF2A
""".strip().split()
promoter_factors = set(x.upper() for x in motif_names)

TSSLoc = namedtuple('TSS', ['chrm', 'strand', 'start', 'stop'])
TSS = namedtuple('TSS', ['loc', 'cnt', 'tpm', 'density'])

VERBOSE = False
QUIET = False
FIX_CHRM_NAMES_FOR_UCSC = False

def log(msg, level=None):
    if QUIET: return
    if level == None or (level == 'VERBOSE' and VERBOSE):
        print >> sys.stderr, msg

def load_TSSs(fp):
    raw_tss_s = []
    total_cnt = 0.0
    for line in fp:
        if line.startswith("track"): continue
        data = line.strip().split("\t")
        loc = TSSLoc(data[0], data[5], int(data[1]), int(data[2]))
        
        cnt = float(data[6])
        total_cnt += cnt
        
        tpm = None
        
        density = numpy.array([float(x) for x in data[-1].split(",")])
        raw_tss_s.append([loc, cnt, tpm, density])
    
    tss_s = []
    for raw_tss in raw_tss_s:
        raw_tss[2] = raw_tss[1]*1e6/total_cnt
        tss_s.append(TSS(*raw_tss))
    
    return tss_s

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Annotate promoters associated with experimental TSSs.')

    parser.add_argument( 'fasta', type=pysam.Fastafile,
        help='Indexed fasta file with the reference sequence.')

    parser.add_argument( 'motifs', type=argparse.FileType("r"),
        help='File containing motifs to analyze.')

    parser.add_argument( 'TSSs', type=argparse.FileType("r"),
        help='Bed file containing annotated TSSs.')
    
    parser.add_argument( '--ucsc', default=False, action='store_true', 
        help='Format the contig names to work with the UCSC genome browser.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', 
                         help='Whether or not to print status information.')
    parser.add_argument( '--quiet', '-q', default=False, action='store_true', 
                         help='Suppress all messages.')

    #parser.add_argument( '--threads', '-t', default=1, type=int,
    #                     help='The number of threads to run.')
    #global NTHREADS
    #NTHREADS = args.threads
        
    args = parser.parse_args()

    global VERBOSE
    if args.verbose: VERBOSE = True 

    global QUIET
    if args.quiet: 
        QUIET = True 
        VERBOSE = False
    
    global FIX_CHRM_NAMES_FOR_UCSC
    FIX_CHRM_NAMES_FOR_UCSC = args.ucsc
    
    return ( args.fasta, args.motifs, args.TSSs)

def main():
    fasta, motifs_fp, tss_s_fp,  = parse_arguments()

    tss_s = load_TSSs(tss_s_fp)
    #for tss in tss_s:
    #    print tss

    motifs = []
    for factor, factor_motifs in load_motifs(
            motifs_fp.name, promoter_factors).iteritems():
        motifs.append(factor_motifs[0])

    return
    """
    for i, tss in enumerate(tss_s):
        if i > 0 and i%1000 == 0: 
            print("Finished {}/{}".format(i, len(tss_s)) )
        seq = tss.get_flanking_seq(fasta)
        for motif, hits in matches.items():
            for hit in hits:
                rv[motif][hit] += 1./len(hits)

    for motif, cnts in rv.items():
        with open("Frac.{}.txt".format(motif), "w") as ofp:
            for cnt in cnts:
                print >> ofp, "{:e}".format(float(cnt)/len(cnts))
    """
    # load the peaks list
    
    # load the fasta file
    
    # load the bam(s)
    pass

if __name__ == '__main__':
    main()
