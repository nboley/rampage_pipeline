#!/usr/bin/env python3
import os, sys, time
import requests, json
from itertools import chain
from collections import defaultdict
import re
from multiprocessing import Value

import numpy

from collections import defaultdict
from itertools import chain

from grit.genes import merge_adjacent_intervals
from grit.files.gtf import load_gtf
from grit.files.reads import RNAseqReads, MergedReads
from grit.frag_len import build_fl_dists_from_annotation
#from grit.f_matrix import build_expected_and_observed_rnaseq_counts, cluster_bins, DesignMatrix, NoObservableTranscriptsError
from grit.f_matrix import build_expected_and_observed_rnaseq_counts, DesignMatrix, NoObservableTranscriptsError
from grit.frequency_estimation import estimate_transcript_frequencies, TooFewReadsError
from grit.transcript import Gene, Transcript
from grit.lib.multiprocessing_utils import ThreadSafeFile, fork_and_wait

# Script to take 1 experiment ID and get to peak calls
# Code adapted from Nathan Boley

# Get fastqs to remap

BASE_URL = "https://www.encodeproject.org"
FASTQ_DIR = "/srv/scratch/dskim89/encodeTSS/fastqs"
BAM_DIR = "/srv/scratch/dskim89/encodeTSS/mapped"
PEAK_DIR = "/srv/scratch/dskim89/encodeTSS/peaks"
WIG_DIR = "/srv/scratch/dskim89/encodeTSS/signal"
GENCODE22 = "/srv/scratch/dskim89/encodeTSS/annotations/gencode.v22.annotation.gtf"
HG22 = "/srv/scratch/dskim89/encodeTSS/annotations/hg22_STAR_2.4.1b_INDEX/"
INCLUDE_LAB_ORIGIN = True

def find_fastqs(experiment_id, is_paired):
    assert isinstance(is_paired, bool)
    URL = "https://www.encodeproject.org/experiments/{}/".format(experiment_id)
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()

    # first build the library mapping
    replicates = {}
    for rep_rec in response_json_dict['replicates']:
        key = (rep_rec['biological_replicate_number'], 
               rep_rec['technical_replicate_number'])
        replicates[key] = rep_rec
    
    if is_paired:
        read_pairs = defaultdict(lambda: [[], []])
    else:
        read_pairs = defaultdict(lambda: [[],])

    for file_rec in response_json_dict['files']:
        file_type = file_rec['file_format']
        if file_type != 'fastq': continue
        rep_key = (file_rec['replicate']['biological_replicate_number'],
                   file_rec['replicate']['technical_replicate_number'] )
        bsid = replicates[rep_key]['library']['biosample']['accession']
        size_range = replicates[rep_key]['library']['size_range']
        depleted_in = ":".join(sorted(
            replicates[rep_key]['library']['depleted_in_term_name']))
        na_term_name = replicates[rep_key]['library']['nucleic_acid_term_name']
        assay = response_json_dict['assay_term_name']
        file_loc = file_rec['href']

        if not is_paired: pair = 0
        else: 
            try: pair = int(file_rec['paired_end'])-1
            except KeyError: return
        
        key = [assay, bsid, size_range, na_term_name, depleted_in]
        if INCLUDE_LAB_ORIGIN: 
            key.append( response_json_dict['lab']['name'] )
        read_pairs[tuple(key)][pair].append(file_loc) # Set up so that first paired end read is in first list, second is in second

        # Also need to extract metadata
    try:
	    control_accession =  response_json_dict['possible_controls'][0]['accession']
    except IndexError:
    	control_accession = None

    return dict(read_pairs), control_accession

def download_fastqs(read_pairs, is_paired, accession, run=False):
	prefixes = []
	bsids = []
	for key in read_pairs.keys():
		[assay, bsid, size_range, na_term_name, depleted_in, author] = key
		prefix = accession + "_" + bsid + "_" + assay + "_"
		prefixes.append(prefix)
		bsids.append(bsid)
		if is_paired:
			for read_pair in range(len(read_pairs[key])):
				for f in read_pairs[key][read_pair]:
					outName = prefix + str(read_pair+1) + ".fastq.gz"
					download = "wget {0}{1} -O {2}/{3}".format(BASE_URL, f, FASTQ_DIR, outName)
					print download
					#if run: os.system(download)
		else:
			pass # All ENCODE data will be paired
	return prefixes, list(set(bsids))

def main():
	# Test with GM12878 (accession: ENCSR000AEI)

	accession = sys.argv[1] # GM12878 RAMPAGE DATA
	is_paired = bool(sys.argv[2])
	threads = sys.argv[3]
	run = True

	# Get FASTQs for RAMPAGE
	print "\nGetting fastqs..."
	[read_pairs_sample, control_accession] = find_fastqs(accession, is_paired)
	[sample_prefixes, sample_bsids] = download_fastqs(read_pairs_sample, is_paired, accession, run)

	[read_pairs_control, placeholder] = find_fastqs(control_accession, is_paired)
	[control_prefixes, control_bsids] = download_fastqs(read_pairs_control, is_paired, accession, run)

	all_prefixes = sample_prefixes + control_prefixes


	# Set up mapping commands 
	print "\nMapping..."
	for prefix in all_prefixes:
		if "RAMPAGE" in prefix:
			if is_paired:
				map_fastq = "bash /users/dskim89/src/map_rampage_star.bash {0} {1}/{2} {3} {4} {5}/{6} {5}/{7}".format(HG22, BAM_DIR, accession, prefix, threads, FASTQ_DIR, prefix + "1.fastq.gz", prefix + "2.fastq.gz")
			else:
				map_fastq = "bash /users/dskim89/src/map_rampage_star.bash {0} {1}/{2} {3} {4} {5}/{6}".format(HG22, BAM_DIR, accession, prefix, threads, FASTQ_DIR, prefix + "1.fastq.gz")
		elif "RNA-seq" in prefix:
			if is_paired:
				map_fastq = "bash /users/dskim89/src/map_rnaseq_star.bash {0} {1}/{2} {3} {4} {5}/{6} {5}/{7}".format(HG22, BAM_DIR, accession, prefix, threads, FASTQ_DIR, prefix + "1.fastq.gz", prefix + "2.fastq.gz")
			else:
				map_fastq = "bash /users/dskim89/src/map_rnaseq_star.bash {0} {1}/{2} {3} {4} {5}/{6}".format(HG22, BAM_DIR, accession, prefix, threads, FASTQ_DIR, prefix + "1.fastq.gz")
		print map_fastq
		#if run: os.system(map_fastq)

	# Index the BAMs
	print "\nIndexing..."
	for bsid in sample_bsids:
		indexRampage = "samtools index {0}/{1}/{1}_{2}_RAMPAGE_markdup.Processed.out.bam {0}/{1}/{1}_{2}_RAMPAGE_markdup.Processed.out.bai".format(BAM_DIR, accession, bsid)
		indexRNA = "samtools index {0}/{1}/{1}_{2}_RNA-seq_Aligned.sortedByCoord.out.bam {0}/{1}/{1}_{2}_RNA-seq_Aligned.sortedByCoord.out.bai".format(BAM_DIR, accession, bsid)
		print indexRampage
		#if run: os.system(indexRampage)
		print indexRNA
		#if run: os.system(indexRNA)

	# Build the wigs. Need to know if stranded or not
	# NOTE: need to add ucsc_tools before building wigs!
	print "\nBuilding wigs..."
	if run: genes = load_gtf(GENCODE22)
	for bsid in sample_bsids:
		print "TESTING MODE: make sure to uncomment grit commands to correctly get strandedness and reverse read strand"
		reads_are_stranded = True
		reverse_read_strand = True

		# First determine strandedness and reverse read strand for RAMPAGE then make bigwig
		filename = "{0}/{1}/{1}_{2}_RAMPAGE_markdup.Processed.out.bam".format(BAM_DIR, accession, bsid)
		if run:
			reads = RNAseqReads(filename).init(ref_genes = genes)
			reads_are_stranded = reads.reads_are_stranded
			reverse_read_strand = reads.reverse_read_strand

		if reads_are_stranded:
			if reverse_read_strand:
				make_bw_rampage = "python /users/dskim89/src/bam2wig.py --mapped-reads-fname {0}/{1}/{1}_{2}_RAMPAGE_markdup.Processed.out.bam -o {3}/{1}_{2}_RAMPAGE --assay rampage -b --ucsc -t {4} -r".format(BAM_DIR, accession, bsid, WIG_DIR, threads)
			else:
				make_bw_rampage = "python /users/dskim89/src/bam2wig.py --mapped-reads-fname {0}/{1}/{1}_{2}_RAMPAGE_markdup.Processed.out.bam -o {3}/{1}_{2}_RAMPAGE --assay rampage -b --ucsc -t {4}".format(BAM_DIR, accession, bsid, WIG_DIR, threads)
		else:
			make_bw_rampage = "python /users/dskim89/src/bam2wig.py --mapped-reads-fname {0}/{1}/{1}_{2}_RAMPAGE_markdup.Processed.out.bam -o {3}/{1}_{2}_RAMPAGE --assay rampage -b --ucsc -t {4} --unstranded".format(BAM_DIR, accession, bsid, WIG_DIR, threads)
		print make_bw_rampage
		if run: os.system(make_bw_rampage)

		# Determine strandedness and reverse read strand for RNA-seq then make bigwig
		filename = "{0}/{1}/{1}_{2}_RNA-seq_Aligned.sortedByCoord.out.bam".format(BAM_DIR, accession, bsid)
		if run:
			reads = RNAseqReads(filename).init(ref_genes = genes)
			reads_are_stranded = reads.reads_are_stranded
			reverse_read_strand = reads.reverse_read_strand		

		if reads_are_stranded:
			if reverse_read_strand:
				make_bw_rnaseq = "python /users/dskim89/src/bam2wig.py --mapped-reads-fname {0}/{1}/{1}_{2}_RNA-seq_Aligned.sortedByCoord.out.bam -o {3}/{1}_{2}_RNA-seq --assay rnaseq -b --ucsc -t {4} -r".format(BAM_DIR, accession, bsid, WIG_DIR, threads)
			else:
				make_bw_rnaseq = "python /users/dskim89/src/bam2wig.py --mapped-reads-fname {0}/{1}/{1}_{2}_RNA-seq_Aligned.sortedByCoord.out.bam -o {3}/{1}_{2}_RNA-seq --assay rnaseq -b --ucsc -t {4}".format(BAM_DIR, accession, bsid, WIG_DIR, threads)
		else:
			make_bw_rnaseq = "python /users/dskim89/src/bam2wig.py --mapped-reads-fname {0}/{1}/{1}_{2}_RNA-seq_Aligned_Aligned.sortedByCoord.out.bam -o {3}/{1}_{2}_RNA-seq --assay rnaseq -b --ucsc -t {4} --unstranded".format(BAM_DIR, accession, bsid, WIG_DIR, threads)
		print make_bw_rnaseq
		if run: os.system(make_bw_rnaseq)

	# Set up peak calling commands. Call peaks for each BSID.
	print "\nCalling peaks..."
	for bsid in sample_bsids:
		rampage_reads = "{0}/{1}/{1}_{2}_RAMPAGE_markdup.Processed.out.bam".format(BAM_DIR, accession, bsid)
		rnaseq_reads = "{0}/{1}/{1}_{2}_RNA-seq_Aligned.sortedByCoord.out.bam".format(BAM_DIR, accession, bsid)
		gff = "{0}/{1}_{2}_grit.gff".format(PEAK_DIR, accession, bsid)
		bed = "{0}/{1}_{2}_peaks.bed".format(PEAK_DIR, accession, bsid)
		call_peaks = "bash /users/dskim89/src/callpeaks_rampage_grit.bash {0} {1} {2} {3} {4} {5}".format(rampage_reads, rnaseq_reads, GENCODE22, gff, bed, threads)
		print call_peaks
		#if run: os.system(call_peaks)

	# Next, use IDR to get reproducible peaks
	print "\nRunning IDR..."
	peakSet1 = "{0}/{1}_{2}_peaks.bed".format(PEAK_DIR, accession, sample_bsids[0])
	peakSet2 = "{0}/{1}_{2}_peaks.bed".format(PEAK_DIR, accession, sample_bsids[1])
	run_idr = "idr --samples {0} {1} --input-file-type bed --rank 7 --output-file {2}/{3}_idr.bed --plot".format(peakSet1, peakSet2, PEAK_DIR, accession)
	print run_idr
	#if run: os.system(run_idr)

	return None

main()