#!/usr/bin/env python3
import os, sys, time

import requests, json

from itertools import chain

from collections import defaultdict

import re

from multiprocessing import Value

MATCH_LAB = True

INDEX_DIR = "/srv/scratch/nboley/RAMPAGE_human/STAR_indexes/STAR_2.4.1b.GRCh38_PhiX_ERCCv2/"
ANNOTATION_GTF = "/mnt/data/annotations/by_release/hg20.GRCh38/GENCODE_ann/gencode.v23/gencode.v23.annotation.gtf"
STAR_CMD = "/usr/local/bin/STAR"

BASE_URL = "https://www.encodeproject.org/"
RNASEQ_STAR_CMD_TEMPLATE = re.sub("\s+", " ", """
%s   --genomeDir %s \
     --readFilesIn {read_fnames}                               \
     --readFilesCommand zcat --runThreadN 24 --genomeLoad LoadAndKeep          \
     --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1\
     --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04             \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000   \
     --outSAMheaderCommentFile COfile.txt                                      \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate                                 \
     --outSAMunmapped Within --outFilterType BySJout                           \
     --outSAMattributes NH HI AS NM MD                                         \
     --outWigType bedGraph --outWigStrand Stranded                             \
     --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000
""" % (STAR_CMD, INDEX_DIR))

RNASEQ_STAR_CMD_TEMPLATE = re.sub("\s+", " ", """
%s   --genomeDir %s \
     --readFilesIn {read_fnames}                               \
     --readFilesCommand zcat --runThreadN 24 --genomeLoad LoadAndKeep          \
     --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1\
     --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04             \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000   \
     --outSAMheaderCommentFile COfile.txt                                      \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate                                 \
     --outSAMunmapped Within --outFilterType BySJout                           \
     --outSAMattributes NH HI AS NM MD                                         \
     --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000
""" % (STAR_CMD, INDEX_DIR)) # SortedByCoordinate --limitBAMsortRAM 30000000000

RAMPAGE_STAR_CMD_TEMPLATE = re.sub("\s+", " ", """
%s    --genomeDir %s
      --readFilesIn {read_fnames}
      --runThreadN 24 
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
      --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000
      --outFilterScoreMinOverLread 0.85 
      --outFilterIntronMotifs RemoveNoncanonicalUnannotated 
      --clip5pNbases 6 15 
      --seedSearchStartLmax 30 """ %  (STAR_CMD, INDEX_DIR))
# --limitBAMsortRAM 10000000000

CAGE_STAR_CMD_TEMPLATE = re.sub("\s+", " ", """
%s    --genomeDir %s
      --readFilesIn {read_fnames}
      --runThreadN 24
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
      --limitBAMsortRAM 60000000000""" % (STAR_CMD, INDEX_DIR))

MARK_PCR_DUP_CMD_TEMPLATE = re.sub("\s+", " ", """
%s    --inputBAMfile {bam_fname}
      --runThreadN 8
      --bamRemoveDuplicatesType UniqueIdentical 
      --runMode inputAlignmentsFromBAM 
      --bamRemoveDuplicatesMate2basesN 15 
      --outFileNamePrefix {op_prefix} 
      --limitBAMsortRAM 60000000000
""" %  STAR_CMD)

GRIT_CMD_TEMPLATE = re.sub("\s+", " ", """
call_peaks.py \
    --rnaseq-reads {RNAseq_reads} --rnaseq-read-type auto
    --rampage-reads {rampage_reads} --rampage-read-type auto
    --reference %s
    --threads 24  
    --exp-filter-fraction 0.10
    --trim-fraction 0.05
    --ucsc --outfname {ofname}
""" % ANNOTATION_GTF)

IDR_CMD_TEMPLATE = re.sub("\s+", " ", """
python3 ~/src/idr/python/idr.py 
    -a {pk1} -b {pk2} --output-file {ofname}
    --peak-merge-method sum --verbose --idr 0.9 --quiet
""")

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
        if MATCH_LAB: 
            key.append( response_json_dict['lab']['name'] )
        read_pairs[tuple(key)][pair].append(file_loc)
        
    return read_pairs

def find_fastqs_from_rampage_experiment(accession):
    return find_fastqs(accession, is_paired=True)

def find_fastqs_from_cage_experiment(accession):
    return find_fastqs(accession, is_paired=False)

def find_matching_rnaseq_experiments(
        bsid, size_range, nucleic_acid, depleted_in, lab=None):
    URL = """https://www.encodeproject.org/search/?type=experiment
             &replicates.library.biosample.accession={bsid}
             &assay_term_name=RNA-seq
             &replicates.library.size_range={size_range}
             &replicates.library.nucleic_acid_term_name={nucleic_acid}"""
    URL = re.sub("\s+", "", URL).format(
        bsid=bsid, size_range=size_range, 
        depleted_in=depleted_in, 
        nucleic_acid=nucleic_acid, 
        lab=lab)
    if depleted_in != '':
        URL += "&replicates.library.depleted_in_term_name={depleted_in}".format(
            depleted_in=depleted_in)
    if lab != None:
        URL += "&lab.name={lab}".format(lab=lab)

    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()['@graph']
    if len(response_json_dict) == 0: 
        print "No results from:", URL
        return None
    try:
        assert len(response_json_dict) == 1
    except: 
        print URL
        raise
    return response_json_dict[0]['accession']

def find_control_experiment(experiment_id):
    URL = "https://www.encodeproject.org/experiments/{}/".format(experiment_id)
    response = requests.get(URL, headers={'accept': 'application/json'})
    pat = '"possible_controls": \["/experiments/(.+?)/"\]'
    controls = sorted(set(re.findall(pat, response.text)))
    assert len(controls) == 1
    return controls[0]

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

    #if os.path.exists(of_prefix): return None, None

    # first make a directory to do everything in, and then downlaod the files
    cmds.append("mkdir {ofname}".format(ofname=of_prefix))
    cmds.append("cd {ofname}".format(ofname=of_prefix))
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
        #cmds.append( MARK_PCR_DUP_CMD_TEMPLATE.format(
        #    bam_fname=initial_bam_fname, op_prefix=bam_fname_prefix ))
        # rampage goes through the PCR duplicate filtering step, so
        # after it finishes remove the old bam, and change the intial
        # bam filename to the merged name
        #cmds.append("rm {}".format(initial_bam_fname))
        #initial_bam_fname = bam_fname_prefix + "Processed.out.bam"
        pass
    
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
            print "FORKED"
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


def build_merged_fname(idr_out_fname):
    experiment_id = idr_out_fname.split(".")[0]
    output_fname = "RAMPAGE_pks.filtered.{}.bed".format(experiment_id)
    with open(idr_out_fname) as fp, open(output_fname, "w") as ofp:
        ofp.write("track type=bed name=RAMPAGE_pks_{}\n".format(experiment_id))
        for i, line in enumerate(fp):
            data = line.split()
            ofp.write("\t".join(
                (data[0], data[1], data[2], 
                 "pk{}".format(i), "1000", data[9])) + "\n")
    return output_fname

def call_peaks_for_experiment(experiment_id):
    # ENCSR000AER
    # ENCSR000AGX
    promoter_fastqs = find_fastqs_from_rampage_experiment(experiment_id)
    #promoter_fastqs = find_fastqs_from_cage_experiment(experiment_id)
    
    mapping_cmds = []
    bam_fnames = {}
    for rampage_key, rampage_fnames in promoter_fastqs.items():
        rampage_bam_fname, rampage_cmd = build_mapping_cmd(
            rampage_key, rampage_fnames)
        if rampage_bam_fname == None: continue
        if rampage_cmd == None: continue
        #exp_id = find_matching_rnaseq_experiments(*rampage_key[1:])
        control_exp_id = find_control_experiment(experiment_id)
                
        for rnaseq_key, rnaseq_fnames in find_fastqs(
                control_exp_id, is_paired=True).items():
            if rampage_key[1:] != rnaseq_key[1:]: continue
            print rampage_key[1:], rnaseq_key[1:]
            rnaseq_bam_fname, rnaseq_cmd = build_mapping_cmd(
                rnaseq_key, rnaseq_fnames)
            if rnaseq_bam_fname != None and rnaseq_cmd != None:
                bam_fnames[rnaseq_key] = rnaseq_bam_fname
                mapping_cmds.append( rnaseq_cmd )
                bam_fnames[rampage_key] = rampage_bam_fname
                mapping_cmds.append( rampage_cmd )

    return mapping_cmds
    #run_in_parallel( mapping_cmds, 16 )
    # call peaks in the matched smaples
    matched_bams = defaultdict(lambda: {})
    for key, fname in bam_fnames.items():
        matched_bams[key[1:]][key[0]] = fname

    peak_calling_cmds = []
    called_peak_fnames = {}
    for key, fnames in matched_bams.items():
        ofname = "RAMPAGE_peaks_{}.narrowPeak".format(
            "_".join(key).replace(">", "GT").replace(" ", ""))
        GRIT_CMD = GRIT_CMD_TEMPLATE.format(
            RNAseq_reads=fnames['RNA-seq'],
            rampage_reads=fnames['RAMPAGE'],
            ofname=ofname)
        peak_calling_cmds.append( GRIT_CMD )
        called_peak_fnames[key] = ofname
 
    #for cmd in peak_calling_cmds:
    #    print cmd
    #run_in_parallel( peak_calling_cmds, 2 )
    """
    # run idr 
    assert len(called_peak_fnames) == 2
    idr_ofname = "{}.idroutput".format(experiment_id)
    cmd = IDR_CMD_TEMPLATE.format(
        pk1=sorted(called_peak_fnames.items())[0][1], 
        pk2=sorted(called_peak_fnames.items())[1][1],
        ofname=idr_ofname)
    #os.system(cmd)
    print cmd
    
    #merged_peaks_fname = build_merged_fname(idr_ofname)
    """
    
    return

def find_finished_rampage_experiments(
        base_dir="/mnt/lab_data/kundaje/projects/TSS_prediction/bams/hg38/"):
    finished_rampage_experiments = set()
    for item in os.listdir(base_dir):
        finished_rampage_experiments.add(item)
    return finished_rampage_experiments

# get all rampage experiments
def find_rampage_experiments():
    URL = "https://www.encodeproject.org/search/?type=experiment&assay_term_name=RAMPAGE&status=released&limit=all"
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()
    biosamples = set()
    for experiment in response_json_dict['@graph']:
        yield experiment['@id'].split("/")[-2]
    return 


if __name__ == '__main__':
    finished_rampage_experiments = find_finished_rampage_experiments()
    mapping_cmds = []
    for exp_id in sorted(set(find_rampage_experiments())):
        if exp_id in finished_rampage_experiments: continue
        res = call_peaks_for_experiment( exp_id )
        if len(res) == 0: continue
        print exp_id, len(res)
        print res[1]
        break
        #mapping_cmds.extend( call_peaks_for_experiment( exp_id ) )
    #run_in_parallel( mapping_cmds, 8 )
    #call_peaks_for_experiment( sys.argv[1] )
