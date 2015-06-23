#!/usr/bin/env python3
import os, sys, time

import requests, json

from itertools import chain

from collections import defaultdict

import re

from multiprocessing import Value

"""
https://www.encodeproject.org/search/?type=experiment&assay_term_name=RNA-seq&limit=all&assembly=hg19&files.file_format=bam&run_type=Paired-ended&replicates.library.size_range=%3E200&replicates.library.nucleic_acid_term_name=polyadenylated%20mRNA
"""

INCLUDE_LAB_ORIGIN = True

INDEX_DIR = "/srv/scratch/nboley/RAMPAGE_human/STAR_index/STARIndex.GRCh38.PhiX.NISC/"
ANNOTATION_GTF = "/data/genomes/GRCh38/gencode.v21.annotation.gtf"
STAR_CMD = "/home/nboley/src/STAR/bin/Linux_x86_64_static/STAR"

BASE_URL = "https://www.encodeproject.org/"
RNASEQ_STAR_CMD_TEMPLATE = re.sub("\s+", " ", """
%s   --genomeDir %s \
     --readFilesIn {read_fnames}                               \
     --readFilesCommand zcat --runThreadN 8 --genomeLoad LoadAndKeep          \
     --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1\
     --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04             \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000   \
     --outSAMheaderCommentFile COfile.txt                                      \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate                                 \
     --outSAMunmapped Within --outFilterType BySJout                           \
     --outSAMattributes NH HI AS NM MD                                         \
     --outWigType bedGraph --outWigStrand Stranded                             \
     --limitBAMsortRAM 30000000000                                             \
     --outSAMtype BAM SortedByCoordinate
""" % (STAR_CMD, INDEX_DIR))

RAMPAGE_STAR_CMD_TEMPLATE = re.sub("\s+", " ", """
%s    --genomeDir %s
      --readFilesIn {read_fnames}
      --runThreadN 8 
      --genomeLoad LoadAndKeep          
      --outSAMunmapped Within 
      --outFilterType BySJout 
      --outSAMattributes NH HI AS NM MD 
      --outFilterMultimapNmax 20 
      --outFilterMismatchNmax 999 
      --outFilterMismatchNoverReadLmax 0.04 
      --alignIntronMin 20 
      --alignIntronMax 1000000 
      --alignMatesGapMax 1000000 
      --alignSJoverhangMin 8 
      --alignSJDBoverhangMin 1 
      --readFilesCommand zcat 
      --outSAMtype BAM SortedByCoordinate 
      --outFilterScoreMinOverLread 0.85 
      --outFilterIntronMotifs RemoveNoncanonicalUnannotated 
      --clip5pNbases 6 15 
      --seedSearchStartLmax 30 
      --limitBAMsortRAM 30000000000""" %  (STAR_CMD, INDEX_DIR))

CAGE_STAR_CMD_TEMPLATE = re.sub("\s+", " ", """
%s    --genomeDir %s
      --readFilesIn {read_fnames}
      --runThreadN 8 
      --genomeLoad LoadAndKeep          
      --outSAMunmapped Within 
      --outFilterType BySJout 
      --outSAMattributes NH HI AS NM MD 
      --outFilterMultimapNmax 20 
      --outFilterMismatchNmax 999 
      --outFilterMismatchNoverReadLmax 0.04 
      --alignIntronMin 20 
      --alignIntronMax 1000000 
      --alignMatesGapMax 1000000 
      --alignSJoverhangMin 8 
      --alignSJDBoverhangMin 1 
      --readFilesCommand zcat 
      --outSAMtype BAM SortedByCoordinate 
      --outFilterScoreMinOverLread 0.85 
      --outFilterIntronMotifs RemoveNoncanonicalUnannotated 
      --clip5pNbases 6 15 
      --seedSearchStartLmax 30 
      --limitBAMsortRAM 30000000000""" % (STAR_CMD, INDEX_DIR))

MARK_PCR_DUP_CMD_TEMPLATE = re.sub("\s+", " ", """
%s    --inputBAMfile {bam_fname}
      --runThreadN 8 
      --bamRemoveDuplicatesType UniqueIdentical 
      --runMode inputAlignmentsFromBAM 
      --bamRemoveDuplicatesMate2basesN 15 
      --outFileNamePrefix {op_prefix} 
      --limitBAMsortRAM 30000000000
""" %  STAR_CMD)

def find_called_peaks(experiment_id, only_merged=True):
    URL = "https://www.encodeproject.org/experiments/{}/".format(experiment_id)
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()

    target = response_json_dict['target']['label']

    # build the replicate mapping
    replicates = {}
    for rep_rec in response_json_dict['replicates']:
        key = (rep_rec['biological_replicate_number'], 
               rep_rec['technical_replicate_number'])
        replicates[key] = rep_rec

    sample_type = response_json_dict['biosample_term_name']
    #sample = replicates.values()[0]['library']['biosample']['aliases']
    #print replicates.values()[0]['library']['biosample']
    #assert False

    for file_rec in response_json_dict['files']:
        file_type = file_rec['file_format']
        if file_type not in ['broadPeak', 'narrowPeak']: continue
        if 'replicate' not in file_rec: 
            rep_key = 'merged'
            bsid = 'merged'
        else:
            if only_merged: continue
            rep_key = (file_rec['replicate']['biological_replicate_number'],
                       file_rec['replicate']['technical_replicate_number'] )
            bsid = replicates[rep_key]['library']['biosample']['accession']
        file_loc = file_rec['href']
        if only_merged:
            yield target, sample_type, file_type, file_loc
        else:
            yield target, sample_type, rep_key, bsid, file_type, file_loc
    return 


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
        read_pairs[tuple(key)][pair].append(file_loc)
        
    return dict(read_pairs)

def find_bams(experiment_id):
    URL = "https://www.encodeproject.org/experiments/{}/".format(experiment_id)
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()
    # first build the library mapping
    replicates = {}
    for rep_rec in response_json_dict['replicates']:
        key = (rep_rec['biological_replicate_number'], 
               rep_rec['technical_replicate_number'])
        replicates[key] = rep_rec
    
    read_pairs = defaultdict(list)
    for file_rec in response_json_dict['files']:
        file_type = file_rec['file_format']
        #if file_type != 'fastq': continue
        if file_type != 'bam': continue
        rep_key = (file_rec['replicate']['biological_replicate_number'],
                   file_rec['replicate']['technical_replicate_number'] )
        bsid = replicates[rep_key]['library']['biosample']['accession']
        size_range = replicates[rep_key]['library']['size_range']
        depleted_in = ":".join(sorted(
            replicates[rep_key]['library']['depleted_in_term_name']))
        na_term_name = replicates[rep_key]['library']['nucleic_acid_term_name']
        assay = response_json_dict['assay_term_name']
        file_loc = file_rec['href']
        key = [assay, bsid, size_range, na_term_name, depleted_in]
        if INCLUDE_LAB_ORIGIN: 
            key.append( response_json_dict['lab']['name'] )
        read_pairs[tuple(key)].append(file_loc)
    
    return dict(read_pairs)

def find_fastqs_from_rampage_experiment(accession):
    return find_fastqs(accession, is_paired=True)

def find_fastqs_from_cage_experiment(accession):
    return find_fastqs(accession, is_paired=False)

def find_fastqs_from_rnaseq_experiment(accession):
    return find_fastqs(accession, is_paired=True)

def build_mapping_cmd(key, fastqs):
    cmds = []
    of_prefix = "_".join(key).replace(">", "GT").replace(" ", "")
    initial_bam_fname = os.path.abspath("./{}/{}".format(
        of_prefix, "Aligned.sortedByCoord.out.bam"))
    bam_fname_prefix = os.path.abspath("./{}/{}".format(
        of_prefix, of_prefix))

    #r1_fastqs, r2_fastqs = [fastqs[0],], [fastqs[1],]

    if key[0] == 'RAMPAGE':
        STAR_CMD_TEMPLATE = RAMPAGE_STAR_CMD_TEMPLATE
    elif key[0] == 'RNA-seq':
        STAR_CMD_TEMPLATE = RNASEQ_STAR_CMD_TEMPLATE
    elif key[0] == 'CAGE':
        STAR_CMD_TEMPLATE = CAGE_STAR_CMD_TEMPLATE
    else: 
        assert False
    
    # first make a directory to do everything in, and then downlaod the files
    cmds.append("mkdir {ofname} && cd {ofname}".format(ofname=of_prefix))
    # now download all of the files, in parallel
    
    for url in chain(*fastqs):
        cmds.append("wget -q {base}{url}".format(base=BASE_URL, url=url))
    # now build the mapping command
    # add the star command, stripping hte raw fname from the full url
    cmds.append( STAR_CMD_TEMPLATE.format(
        read_fnames=" ".join((",".join(x.split("/")[-1] for x in r_fastqs) 
                              for r_fastqs in fastqs ))
    ))
    
    if key[0] == 'RAMPAGE':
        cmds.append( MARK_PCR_DUP_CMD_TEMPLATE.format(
            bam_fname=initial_bam_fname, op_prefix=bam_fname_prefix ))
        # rampage goes through the PCR duplicate filtering step, so
        # after it finishes remove the old bam, and change the intial
        # bam filename to the merged name
        cmds.append("rm {}".format(initial_bam_fname))
        initial_bam_fname = bam_fname_prefix + "Processed.out.bam"

    cmds.append("mv {} {}.bam".format(initial_bam_fname, bam_fname_prefix))
    # add the index cmd
    cmds.append("sambamba index -t 16 {}.bam".format(bam_fname_prefix))
    
    return bam_fname_prefix + ".bam", " &&\n".join(cmds)

def run_in_parallel(cmds_queue, NPROC):
    running_ps = Value("i")
    while len(cmds_queue) > 0:
        # if we already have the maximum number of processes running,
        # then sleep and wait for a process to become free
        if running_ps.value >= NPROC:
            time.sleep(1)
            continue
        
        # get commands to run, and increment the number of processes
        # if we can't then we are out of commands, so break
        try: cmds = cmds_queue.pop()
        except KeyError: break
        
        with running_ps.get_lock(): 
            running_ps.value += 1
        # fork a process. If we are still in main, then do nothing. Otherwise,
        # run the grabbed commands, decrement the running processes value
        # and exit
        pid = os.fork()
        if pid != 0:
            continue
        else:
            print(cmds)
            os.system(cmds)
            with running_ps.get_lock(): 
                running_ps.value -= 1
            os._exit(0)
    
    # wait for outstanding processes to finish
    while True:
        with running_ps.get_lock():
            if running_ps.value == 0: break
        time.sleep(1)
        continue

    return

def find_rnaseq_experiments():
    URL = "https://www.encodeproject.org/search/?type=experiment&assay_term_name=RNA-seq&limit=all&assembly=hg19&files.file_format=bam&run_type=Paired-ended&replicates.library.size_range=%3E200&replicates.library.nucleic_acid_term_name=polyadenylated%20mRNA"
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()
    biosamples = set()
    for experiment in response_json_dict['@graph']:
        yield experiment['@id'].split("/")[-2]
    return 

def find_chipseq_experiments():
    URL = "https://www.encodeproject.org/search/?type=experiment&assay_term_name=ChIP-seq&assembly=mm9&target.investigated_as=transcription%20factor&limit=all&format=json"
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()
    biosamples = set()
    for experiment in response_json_dict['@graph']:
        yield experiment['@id'].split("/")[-2]
    return 

def find_rnaseq_fastqs_and_bams():
    for exp_id in find_rnaseq_experiments():
        fastqs = find_fastqs_from_rnaseq_experiment(exp_id)
        if fastqs == None: continue
        print "="*80
        print "fastqs"
        print fastqs
        print "bams"
        print find_bams(exp_id )

def build_tf_list():
    with open("mouse_tfs.txt", "w") as ofp:
        chipseq_exps = list(find_chipseq_experiments())
        for i, exp_id in enumerate(chipseq_exps):
            for res in find_called_peaks(exp_id, True):
                target, sample_type, file_type, file_loc = res
                human_readable_ofname = (
                    "_".join((target, sample_type.replace(" ", "-")))
                    + ".%s.bed" % file_type )
                ofname = file_loc.split("/")[-1]
                print "Downloading %i/%i" % (i+1, len(chipseq_exps))
                cmd = "wget --quiet {URL}; bigBedToBed {FNAME} {HR_FNAME}".format(
                    URL=BASE_URL+file_loc, 
                    FNAME=ofname, HR_FNAME=human_readable_ofname )
                os.system(cmd)
                print >> ofp, "\t".join(
                    (exp_id, target, sample_type, 
                     BASE_URL[:-2]+file_loc, 
                     human_readable_ofname))
    
    return
    
if __name__ == '__main__':
    call_peaks_for_experiment( sys.argv[1] )
