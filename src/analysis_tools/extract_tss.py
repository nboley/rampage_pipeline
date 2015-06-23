import os, sys

try: import grit
except ImportError: sys.path.insert(0, "/home/nboley/grit/grit/")

from grit.files.gtf import load_gtf

genes = load_gtf(sys.argv[1])
for gene in genes:
    promoters = gene.extract_elements()['promoter']
    for start, stop in promoters:
        print "\t".join((gene.chrm, str(start), str(stop), 
                        "tss", "1000", gene.strand))
