__version__ = 0.1

import os, sys
import subprocess

from collections import defaultdict, namedtuple

TSSMatch = namedtuple('TSSMatch', ['score', 'dist', 'pos', 'ref_pos'])

import numpy
from scipy.signal import convolve

import pybedtools

from grit.files.gtf import load_gtf
from grit.files.reads import fix_chrm_name_for_ucsc

import gzip

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_annotation_tss_s(fname):
    if fname.endswith('.gz'):
        fp = gzip.open(fname)
    else:
        fp = open(fname)
    
    genes = load_gtf(fp)
    tss_s = []
    for gene in genes:
        for start, stop in gene.extract_elements()['promoter']:
            tss_s.append(
                (gene.chrm, str(start), str(stop), gene.strand))
    return tss_s

def write_tss_s_to_bed(tss_s, ofname, element_name='TSS'):
    with open(ofname, "w") as ofp: 
        for chrm, start, stop, strand in tss_s:
            ofp.write("{}\t{}\t{}\t{}\t1000\t{}\n".format(
                fix_chrm_name_for_ucsc(chrm), 
                start, stop, element_name, strand))
    return

def compare_promoters(ref_fname, other_fname, filetype=None, min_score=None):
    cmd ="bedtools closest -b {} -a {} -s -d -t all".format(
        ref_fname, other_fname)
    res = subprocess.check_output(cmd, shell=True)
    
    # find the index of the start of the gencode contig, and try
    # and geus the filetype if it's not set
    second_column_start = None
    for line in res.split("\n"):
        # skiup header lines
        if line.startswith('track'): continue
        
        data = line.split("\t")
        # if we don't have a filetype, try and infer it from the line length
        if filetype == None:
            if len(data) == 18: 
                filetype = 'bed'
            else:
                filetype = 'idr'
        # try to find the index of the contig. If we can't find one, 
        # try the next line
        try: 
            second_column_start = data[1:].index(data[0])+1
        except ValueError, inst:
            continue
        break
    assert second_column_start != None

    if filetype == 'idr' and min_score == None:
        min_score = 1.30103 # -log10(0.05)
    
    scores_and_dist = {}
    for line in res.split("\n"):
        if line.startswith('track'): continue
        data = line.strip().split("\t")
        if data == ['']: continue
        if data[second_column_start] == '.': continue
        dist = int(data[-1])
        # if this is raw grit output
        if filetype == 'bed':
            score = float(data[6])
        # if this is IDR output
        elif filetype == 'idr':
            # mean signal score
            cnt_score = float(data[10]) + float(data[13])
            # IDR score
            idr_score = float(data[7])
            score = (idr_score, cnt_score)
        else:
            assert False
            
        # skip lines with too low of a score
        if min_score != None and score[0] < min_score:
            continue
        
        # when we detect a duplicate, always match the highest score
        key = tuple(data[second_column_start:second_column_start+3])
        if key not in scores_and_dist: 
            scores_and_dist[key] = TSSMatch(
                score, dist, 
                "_".join(map(str, data[:3])), 
                "_".join(map(str, data[second_column_start:second_column_start+3])))
        elif scores_and_dist[key][0] < score:
            scores_and_dist[key] = TSSMatch(
                score, dist, 
                "_".join(map(str, data[:3])), 
                "_".join(map(str, data[second_column_start:second_column_start+3])))

    
    for x in sorted(scores_and_dist.values(), 
                    key=lambda x:x.score, reverse=True):
        yield x

def load_annotation(fname):
    ref_tss_fname = "TSS." + os.path.basename(fname).replace(".gz", "")
    try: 
        with open(ref_tss_fname) as fp:
            return ref_tss_fname
    except:
        pass
    tss_s = load_annotation_tss_s(fname)
    write_tss_s_to_bed(tss_s, ref_tss_fname, "GENCODE_TSS")
    return ref_tss_fname

def build_TSS_distance_plot(
        tss_matches, fname, smooth_window_size=100, max_distance=200):
    fig = plt.figure( num=None, figsize=(9, 9))
    plt.title("Score vs Frac of TSSs within %i BP of GENCODE TSS" 
              % max_distance)
    plt.xlabel("Peak Rank")
    plt.ylabel("Frac of TSSs within %i BP of GENCODE TSS" % max_distance)
    plt.ylim((0,1))
    scores = numpy.array([int(x.dist < max_distance) for x in tss_matches])
    window = numpy.ones(smooth_window_size, dtype=float)/smooth_window_size
    smooth_scores = convolve(scores, window, mode='valid')
    plt.plot(smooth_scores)
    plt.savefig(fname)
    return

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
Program: TSS Peak Calling QC Code
Version: {PACKAGE_VERSION}
Contact: Nathan Boley <npboley@gmail.com>
""".format(PACKAGE_VERSION=__version__))

    parser.add_argument( '--annotation', type=file, required=True,
                         help='Annotation in GTF format.')
    parser.add_argument( '--GRIT-bed-peaks', nargs="*", type=file, default=[],
        help='Bed file containg TSS peaks called by GRIT.')
    parser.add_argument( '--IDR-bed-peaks', nargs="?", type=file,
        help='IDR bed output file.')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
        
    gencode_tss_fname = load_annotation(args.annotation.name)

    for raw_peaks_fp in args.GRIT_bed_peaks:
        tss_matches =  list(compare_promoters(
            gencode_tss_fname, raw_peaks_fp.name))
        build_TSS_distance_plot(
            tss_matches, 
            os.path.basename(raw_peaks_fp.name)+".png", 
            200, 100)

    if args.IDR_bed_peaks != None:
        tss_matches =  list(compare_promoters(
            gencode_tss_fname, args.IDR_bed_peaks.name))
        build_TSS_distance_plot(
            tss_matches, 
            os.path.basename(args.IDR_bed_peaks.name)+".png", 
            200, 100)

if __name__ == '__main__':
    main()
