import os, sys

from collections import namedtuple, OrderedDict, defaultdict

import re

import numpy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pysam

sys.path.insert(0, "/users/nboley/src/TF_binding/")
#from motif_tools import Motif, iter_motifs
from load_motifs import load_motifs 

USE_START = False
USE_MAX = True
USE_WEIGHTED_AVERAGE = False

VERBOSE = False
QUIET = False
FIX_CHRM_NAMES_FOR_UCSC = False

motif_names = """
TATA
BDP1
HMGN3
GTF2I
GTF2A
""".strip().split()
promoter_factors = set(x.upper() for x in motif_names)

try: 
    rev_comp_table = str.maketrans("ACGT", "TGCA")
except:
    import string
    rev_comp_table = string.maketrans("ACGT", "TGCA")

def log(msg, level=None):
    if QUIET: return
    if level == None or (level == 'VERBOSE' and VERBOSE):
        print >> sys.stderr, msg

TSSLoc = namedtuple('TSS', ['chrm', 'strand', 'start', 'stop'])
TSSData = namedtuple('TSS', ['loc', 'cnt', 'tpm', 'density'])
class TSS(TSSData):
    def get_sequence(self, fasta, flank_size=0):
        region = (self.loc.chrm, 
                  self.loc.start-flank_size, 
                  self.loc.stop+flank_size)
        seq = str(fasta.fetch(*region).upper())
        if self.loc.strand == '-':
            return seq.translate(rev_comp_table)[::-1] 
        else:
            return seq
    
    def __len__(self):
        return self.loc.stop - self.loc.start

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

def build_plot_for_factor(tss_s, motif, fasta, flank_length):
    fasta = pysam.Fastafile(fasta.filename)
    agg_score = numpy.zeros(flank_length*2-len(motif), dtype=float)
    for tss_i, tss in enumerate(tss_s):
        if tss_i%1000 == 0: print motif.name, tss_i, "/", len(tss_s)
        try: seq = tss.get_sequence(fasta, flank_length)
        except IndexError: continue
        if seq == '': continue
        scores = numpy.array([
            score for pos, score in motif.iter_pwm_score(seq)])
        if USE_WEIGHTED_AVERAGE:
            for offset, pos_cnt in zip(xrange(len(tss)), tss.density):
                agg_score[:] += (float(pos_cnt)/tss.cnt
                             )*scores[offset:offset+2*flank_length-len(motif)]
        elif USE_MAX:
            offset = numpy.argmax(tss.density)
            agg_score[:] += scores[offset:offset+2*flank_length-len(motif)]
        elif USE_START:
            offset = 0
            agg_score[:] += scores[offset:offset+2*flank_length-len(motif)]
        else:
            assert False
        #print (scores.sum()/len(motif))/len(scores), len(scores)-len(tss)+len(motif)
        #print seq
    
    plt.plot(len(agg_score)*agg_score/agg_score.sum())
    plt.savefig(motif.name + ".png")
    return

def main():
    fasta, motifs_fp, tss_s_fp,  = parse_arguments()

    motifs = []
    for factor, factor_motifs in load_motifs(
            motifs_fp.name, promoter_factors).iteritems():
        motifs.append(factor_motifs[0])

    tss_s = load_TSSs(tss_s_fp)

    # the amount of flanking sequence we will grab, on each side 
    flank_length = 50
    pids = []
    for motif in motifs:
        pid = os.fork()
        #pid = 0
        if pid == 0:
            build_plot_for_factor(tss_s, motif, fasta, flank_length)
            os._exit(0)
        else:
            pids.append(pid)
    
    for pid in pids:
        os.waitpid(pid, 0)

    return

if __name__ == '__main__':
    main()
