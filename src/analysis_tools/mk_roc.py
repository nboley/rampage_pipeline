import os, sys
import subprocess

from collections import defaultdict

o = 0 #4

def main():
    fname = sys.argv[1]
    cmd ="bedtools closest -b gencode.v19.TSS.notlow.chr1.gff -a {} -t first".format(
        fname)
    cmd ="bedtools closest -b gencode.v19.TSS.notlow.gff -a {} -t first".format(
        fname)
    res = subprocess.check_output(cmd, shell=True)
    scores_and_dist = {}
    for line in res.split("\n"):
        data = line.split()
        if data == []: continue
        if data[10] == '.': continue
        dist = min(abs(int(data[1])-int(data[13-o])), 
                   abs(int(data[2])-int(data[13-o])),
                   abs(int(data[1])-int(data[14-o])),
                   abs(int(data[2])-int(data[14-o])),
        )
        score = float(data[6]) #6 for GRIT #4 for bed score
        key = tuple(data[10-o:15-o])
        if key not in scores_and_dist: 
            scores_and_dist[key] = (score, dist, "_".join(map(str, data[:3])), "_".join(map(str, data[10-o:15-o])))
        elif scores_and_dist[key][0] < score:
            scores_and_dist[key] = (score, dist, "_".join(map(str, data[:3])), "_".join(map(str, data[10-o:15-o])))
    for x in sorted(scores_and_dist.values(), key=lambda x:x[0], reverse=True):
        print "\t".join(map(str,x))
        
main()
