#!/usr/bin/env python3
import os, sys, time

import requests, json

import sqlite3

from itertools import chain

import re

from collections import namedtuple, defaultdict
def namedtuple_factory(cursor, row):
    """
    Usage:
    con.row_factory = namedtuple_factory
    """
    fields = [col[0] for col in cursor.description]
    Row = namedtuple("Row", fields)
    return Row(*row)

from multiprocessing import Value

BASE_URL = "https://www.encodeproject.org/"
STAR_CMD_TEMPLATE = re.sub("\s+", " ", """
STAR --genomeDir /srv/scratch/nboley/RAMPAGE_human/STAR_index/STARIndex.hg19.PhiX.NISC/ \
     --readFilesIn {read1_fnames} {read2_fnames}                               \
     --readFilesCommand zcat --runThreadN 16 --genomeLoad LoadAndKeep          \
     --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1\
     --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04             \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000   \
     --outSAMheaderCommentFile COfile.txt                                      \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate                                 \
     --outSAMunmapped Within --outFilterType BySJout                           \
     --outSAMattributes NH HI AS NM MD                                         \
     --outWigType bedGraph --outWigStrand Stranded                             \
     --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM
""")

GRIT_CMD_TEMPLATE = re.sub("\s+", " ", """
python ~/grit//grit/bin/call_peaks.py 
    --reference /data/genomes/hg19/gencode.v19.annotation.gtf 
    --rnaseq-reads {RNAseq_reads} 
    --rampage-reads {rampage_reads} 
    --threads 8  
    --exp-filter-fraction 0.10
    --trim-fraction 0.05
    --ucsc --out-fname {ofname}
""")

IDR_CMD_TEMPLATE = re.sub("\s+", " ", """
python3 ~/src/idr/python/idr.py 
    -a {pk1} -b {pk2} --output-file {ofname}
    --peak-merge-method sum --verbose --idr 0.10 --quiet
""")

FastqRecord = namedtuple('FastqRecord', [
    'bsid', 'assay_type', 'biorep', 'read_pair', 'url' ])

# get all rampage experiments
def find_rampage_biosamples():
    URL = "https://www.encodeproject.org/search/?type=experiment&assay_term_name=RAMPAGE&status=released&files.file_format=bam&frame=embedded"
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()
    biosamples = set()
    for experiment in response_json_dict['@graph']:
        for rep in experiment['replicates']:
            #print json.dumps(rep, indent=4, separators=(',', ': '))
            biosamples.add( rep['library']['biosample']["accession"] )
    return biosamples

def find_fastqs_from_biosample(bsid):
    URL = "https://www.encodeproject.org/search/?searchTerm={}&type=experiment&frame=embedded".format(bsid)
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()
    for experiment in response_json_dict['@graph']:
        for file_rec in experiment['files']:
            # skiup files that aren't fastqs
            file_type = file_rec['file_format']
            if file_type != 'fastq': continue# not in ('fastq', 'bam'): continue
            # skip replicates that aren't paired
            try: rep = file_rec['replicate']
            except: continue
            
            if not rep['paired_ended']: continue

            yield FastqRecord(*[bsid, 
                                experiment['assay_term_name'], 
                                rep['biological_replicate_number'], 
                                int(file_rec['paired_end']), 
                                file_rec['href']])
    return 

class SamplesDB:
    def _init_fastq_files_table(self):
        c = self.conn.cursor()
        try: 
            if self.rebuild:
                c.execute('DROP TABLE fastq_files;');
            c.execute('''
        CREATE TABLE fastq_files
        (    
            bsid text, 
            assay_type text, 
            biorep int, 
            read_pair int,
            url text PRIMARY KEY,
            location text,
            mapped bool
        );
            ''')
        # if the table already exists, don't create it
        except sqlite3.OperationalError:
            if self.rebuild: raise
            else: self.conn.rollback()
        else:
            self.conn.commit()
        c.close()
        return

    def _init_bams_table(self):
        c = self.conn.cursor()
        try: 
            if self.rebuild:
                c.execute('DROP TABLE bam_files;');
            c.execute('''
        CREATE TABLE bam_files
        (    
            bsid text, 
            biorep int, 
            assay_type text, 
            location text PRIMARY KEY
        );
            ''')
        # if the table already exists, don't create it
        except sqlite3.OperationalError:
            if self.rebuild: raise
            else: self.conn.rollback()
        else:
            self.conn.commit()
        c.close()
        return
    
    def __init__(self, fname, rebuild=False):
        self.conn = sqlite3.connect(fname)
        self.rebuild = rebuild
        self.conn.row_factory = namedtuple_factory
        self._init_fastq_files_table()
        self._init_bams_table()
    
    def add_fastq(self, values):
        with self.conn:
            c = self.conn.cursor()
            c.execute('INSERT INTO fastq_files VALUES (?, ?, ?, ?, ?, 0, 0);', 
                      values);
        return

    def add_bam(self, bsid, biorep, assay_type, fname):
        with self.conn:
            c = self.conn.cursor()
            c.execute('INSERT INTO bam_files VALUES (?, ?, ?, ?);', 
                      (bsid, biorep, assay_type, fname));
        return

    def get_grpd_fastqs(self):
        bsid_biorep = defaultdict(lambda: {1: [], 2: []})
        with self.conn:
            c = self.conn.cursor()
            # hope that the paired fastqs are next to each other in ordered 
            # urls - I don't see an easier way
            for r in c.execute('''SELECT * FROM fastq_files 
                                  ORDER BY bsid, biorep, assay_type, url'''):
                key = (r.bsid, r.biorep, r.assay_type)
                bsid_biorep[key][r.read_pair].append(r.url)
        
        return dict(bsid_biorep)

    def get_grpd_bams(self):
        bsid_biorep = defaultdict(lambda: {'RNA-seq': [], 'RAMPAGE': []})
        with self.conn:
            c = self.conn.cursor()
            # hope that the paired fastqs are next to each other in ordered 
            # urls - I don't see an easier way
            for r in c.execute('''SELECT * FROM bam_files 
                                  ORDER BY bsid, biorep, assay_type, location'''):
                bsid_biorep[(r.bsid, r.biorep)][r.assay_type].append(r.location)
        
        return dict(bsid_biorep)

def build_mapping_cmd(bsid, biorep, assay_type, fastqs):
    cmds = []
    of_prefix = "{}_{}_{}".format(bsid, biorep, assay_type)
    bam_fname = os.path.abspath("./{}/{}".format(
        of_prefix, "Aligned.sortedByCoord.out.bam"))
    r1_fastqs, r2_fastqs = fastqs[1], fastqs[2]

    # first make a directory to do everything in, and then downlaod the files
    cmds.append("mkdir {ofname} && cd {ofname}".format(ofname=of_prefix))
    # now download all of the files, in parallel
    for url in chain(r1_fastqs, r2_fastqs):
        cmds.append("wget -q {base}{url}".format(base=BASE_URL, url=url))
    # now build the mapping command
    # add the star command, stripping hte raw fname from the full url
    cmds.append( STAR_CMD_TEMPLATE.format(
        read1_fnames=",".join(x.split("/")[-1] for x in r1_fastqs), 
        read2_fnames=",".join(x.split("/")[-1] for x in r2_fastqs)) )
                              
    # add the index cmd
    cmds.append("sambamba index -t 16 Aligned.sortedByCoord.out.bam")
    
    return bam_fname, " &&\n".join(cmds)

def add_fastqs_to_db(db):
    for bsid in find_rampage_biosamples():
        for record in find_fastqs_from_biosample(bsid):
            try: db.add_fastq( record )
            except sqlite3.IntegrityError:
                print("Record already exists:", record)
    return

def add_bams_to_db(db, bam_fnames):
    for (bsid, biorep, assay_type), fname in bam_fnames.items():
        try: 
            db.add_bam(bsid, biorep, assay_type, fname)
        except sqlite3.IntegrityError:
            print("BAM already exists:", fname)
    
    return

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

def download_and_map(db, NPROC=4):
    cmds_queue = []
    bam_fnames = {}
    for (bsid, biorep, assay_type), fastqs in db.get_grpd_fastqs().items():
        bam_fname, cmds = build_mapping_cmd(bsid, biorep, assay_type, fastqs)
        cmds_queue.append(cmds)
        bam_fnames[(bsid, biorep, assay_type)] = bam_fname
    
    run_in_parallel( cmds_queue, NPROC )
    
    return bam_fnames

def call_peaks(db, NPROC=8):
    cmds_queue = []
    peak_fnames = defaultdict(dict)
    bsid_merged_fnames = defaultdict(
        lambda: {'num_reps': 0, 'RNA-seq': [], 'RAMPAGE': []})
    for (bsid, biorep), fnames in db.get_grpd_bams().items():
        RNAseq_fnames = fnames['RNA-seq']
        RAMPAGE_fnames = fnames['RAMPAGE']
        if len(RNAseq_fnames) == 0 or len(RAMPAGE_fnames) == 0: 
            continue
        assert len(RNAseq_fnames) == 1 and len(RAMPAGE_fnames) == 1
        
        bsid_merged_fnames[bsid]['num_reps'] += 1
        bsid_merged_fnames[bsid]['RNA-seq'].append(RNAseq_fnames[0])
        bsid_merged_fnames[bsid]['RAMPAGE'].append(RAMPAGE_fnames[0])
        
        ofname = "RAMPAGE_peaks_{}_{}.narrowPeak".format(bsid, biorep)
        GRIT_CMD = GRIT_CMD_TEMPLATE.format(
            RNAseq_reads=RNAseq_fnames[0],
            rampage_reads=RAMPAGE_fnames[0],
            ofname=ofname)
        cmds_queue.append(GRIT_CMD)
        peak_fnames[bsid][biorep] = ofname
    
    run_in_parallel( cmds_queue, NPROC )
    #for key, vals in bsid_merged_fnames.items():
    #    print( key, vals )
    
    return peak_fnames

def run_idr(peak_fnames, NPROC=32):
    cmds_queue = []
    ofnames = {}
    for bsid, fnames in peak_fnames.items():
        if len(fnames) == 1: continue
        assert set(fnames.keys()) == set((1,2)), fnames
        ofname = "{}.idroutput".format(bsid)
        cmd = IDR_CMD_TEMPLATE.format(
            pk1=fnames[1], pk2=fnames[2], ofname=ofname)
        cmds_queue.append(cmd)
        ofnames[ofname] = (fnames[1], fnames[2])
        
    run_in_parallel( cmds_queue, NPROC )
    return ofnames

def build_merged_fname(idr_out_fname, input_fnames):
    bsid = idr_out_fname.split(".")[0]
    output_fname = "RAMPAGE_pks.filtered.{}.bed".format(bsid)
    with open(idr_out_fname) as fp, open(output_fname, "w") as ofp:
        ofp.write("track type=bed name=RAMPAGE_pks_{}\n".format(bsid))
        for i, line in enumerate(fp):
            data = line.split()
            ofp.write("\t".join(
                (data[0], data[1], data[2], 
                 "pk{}".format(i), "1000", data[9])) + "\n")
    return output_fname

def build_merged_fnames(idr_fnames):
    for idr_out_fname, input_fnames in idr_fnames.items():
        try: build_merged_fname(idr_out_fname, input_fnames)
        except: print( "ERROR: Can't build {}".format(idr_out_fname) )
    return

def main():
    # load the sql lite db
    db = SamplesDB('samples.db', False)
    #populate_db(db)
    #add_bams_to_db(db, download_and_map(db))
    peak_fnames = call_peaks(db)
    idr_fnames = run_idr(peak_fnames)
    build_merged_fnames( idr_fnames )
    return

if __name__ == '__main__':
    main()
