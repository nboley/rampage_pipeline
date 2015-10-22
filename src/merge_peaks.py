import os, sys
import cPickle as pickle
from collections import defaultdict
import re
import urllib

import numpy

from grit.files.gtf import flatten

PEAKS_PATH = "/mnt/lab_data/kundaje/projects/TSS_prediction/encodeTSS/peaks/"
PEAKS_PATH = "/home/nboley/Desktop/rampage_paper/peaks/"

trim_fraction = 0.05

def iter_bed_peaks(fname):
    with open(fname) as fp:
        for line in fp:
            if line.startswith('track'): continue
            data = line.split()
            contig = data[0]
            start = int(data[1])
            stop = int(data[2])
            strand = data[5]
            yield (contig, strand, start, stop), data
    return

def load_idr_peaks(path=PEAKS_PATH):
    peaks = defaultdict(list)
    for i, fname in enumerate(os.listdir(path)):
        if not fname.endswith("_idr.bed"): continue
        experiment = fname.split("_")[0]
        assay = get_assay_name(experiment)
        print >> sys.stderr, assay, fname
        if assay == 'CAGE': continue
        for (contig, strand, start, stop), data in iter_bed_peaks(
                os.path.join(path, fname)):
            peaks[(contig, strand)].append((start, stop))
    return peaks

def load_sample_peaks(path=PEAKS_PATH):
    counts = defaultdict(float)
    peaks = defaultdict(list)
    for fname in os.listdir(path):
        print fname
        if not fname.endswith("_peaks.bed"): continue
        for (contig, strand, start, stop), data in iter_bed_peaks(
                os.path.join(path, fname)):
            peaks[(contig, strand)].append((start, stop, data, fname))
            counts[fname] += float(data[6])
        #if len(counts) > 15: break
    return peaks, dict(counts)

def group_peaks_in_contig(contig_peaks):
    grpd_intervals = [[],]
    curr_start, curr_stop = contig_peaks[0][:2]
    for start, stop, data, fname in contig_peaks:
        if start < curr_stop:
            curr_stop = max(stop, curr_stop)
            grpd_intervals[-1].append(((start, stop), data, fname))
        else:
            curr_start, curr_stop = start, stop
            grpd_intervals.append([((start, stop), data, fname),])
    return grpd_intervals

def merge_peaks(peaks):
    merged_peaks = {}
    for key, contig_peaks in peaks.iteritems():
        intervals = [x[:2] for x in contig_peaks]
        merged_peaks[key] = [(x[0], x[1], None, 'IDR')
                             for x in flatten(intervals)]
    return merged_peaks

def trim_coverage(cov_region, start, stop):
    total_cov = cov_region.sum()
    cov_cumsum = cov_region.cumsum()-cov_region[0]
    try: trim_start = numpy.flatnonzero(
            cov_cumsum < int(trim_fraction*total_cov)).max()
    except:
        trim_start = 0
    try: trim_stop = numpy.flatnonzero(
            cov_cumsum > (1.0-trim_fraction)*total_cov).min()
    except: trim_stop=len(cov_region)-1
    while trim_start < len(cov_region)-1 and cov_region[trim_start] == 0:
        trim_start += 1
    while trim_stop > trim_start and cov_region[trim_stop] == 0:
        trim_stop -= 1
    return trim_start+start, trim_stop+start, cov_region[trim_start:trim_stop+1]


def extract_coverage(data):
    return numpy.array([float(x) for x in data[-1].split(",")], dtype=float)

def iter_merged_peaks(grpd_peaks, depths):
    for peaks in grpd_peaks:
        idr_peaks = [pk for pk in peaks if pk[2] == 'IDR']
        for idr_peak in idr_peaks:
            pk_start, pk_stop = idr_peak[0]

            # build the read coverage
            cov = numpy.zeros(pk_stop - pk_start + 1, dtype=float)
            for peak in peaks:
                if peak[2] == 'IDR': continue
                if pk_start > peak[0][0] or pk_stop < peak[0][1]:
                    continue
                scale = 1e6/(depths[peak[-1]]*max(5, len(peaks)))
                cov[peak[0][0]-pk_start:peak[0][1]-pk_start] += (
                    scale*extract_coverage(peak[1]) )
            
            pk_start, pk_stop, cov = trim_coverage(cov, pk_start, pk_stop)
            yield pk_start, pk_stop, cov
    return

def print_merged_peaks(merged_peaks):
    for (contig, strand), intervals in merged_peaks.iteritems():
        if contig.startswith('chrGL'): contig = contig[3:]
        if contig.startswith('chrKI'): contig = contig[3:]
        if contig.startswith('chrKN'): contig = contig[3:]
        if contig.startswith('chrJH'): contig = contig[3:]
        if contig == 'chrphiX174': continue
        for start, stop in intervals:
            print "%s\t%i\t%i\t.\t1000\t%s" % (
                contig, start, stop, strand)
    return

def clean_contig_name(contig):
    if contig.startswith('chrGL'): contig = contig[3:]
    if contig.startswith('chrKI'): contig = contig[3:]
    if contig.startswith('chrKN'): contig = contig[3:]
    if contig.startswith('chrJH'): contig = contig[3:]
    return contig

def get_assay_name(experiment):
    cage_pat = 'assay_term_name": "CAGE"'
    rampage_pat = 'assay_term_name": "RAMPAGE"'
    url = "https://www.encodeproject.org/experiments/%s/?format=json" % experiment
    data = urllib.urlopen(url).read()
    cage_res = re.findall(cage_pat, data)
    rampage_res = re.findall(rampage_pat, data)
    assert len(cage_res) == 0 or len(rampage_res) == 0
    if len(cage_res) > 0:
        return 'CAGE'
    if len(rampage_res) > 0:
        return 'RAMPAGE'
    assert False

def main():
    with open("idr_peaks.rampage.obj") as fp:
        idr_peaks = pickle.load(fp)
    #idr_peaks = merge_peaks(load_idr_peaks())
    #with open("idr_peaks.rampage.obj", "w") as ofp:
    #    pickle.dump(idr_peaks, ofp)
    print >> sys.stderr, "Finished Loading IDR Peaks"

    with open("sample_peaks.obj") as fp:
        sample_peaks, depths = pickle.load(fp)
    #sample_peaks, depths = load_sample_peaks()
    #with open("sample_peaks.obj", "w") as ofp:
    #    pickle.dump((sample_peaks, depths), ofp)
    print >> sys.stderr, "Finished Loading Sample Peaks"

    for contig, strand in idr_peaks.keys():
        print >> sys.stderr, contig, strand
        if contig == 'chrphiX174': continue
        if contig != clean_contig_name(contig): continue
        grpd_peaks = group_peaks_in_contig(
            sorted(idr_peaks[(contig, strand)] + sample_peaks[(contig, strand)])
        )
        for start, stop, cov in iter_merged_peaks(grpd_peaks, depths):
            score = min(1000, int(cov.sum()*100))
            print "%s\t%i\t%i\tNone\t%i\t%s" % (
                clean_contig_name(contig), start, stop, score, strand)
    
    return

if __name__ == '__main__':
    main()
