import os, sys
with open(sys.argv[1]) as fp:
    for i, line in enumerate(fp):
        if line.startswith("track"): 
            print line,
            continue
        data = line.split()
        new_data = data[:4] + [str(int(float(data[6]))), data[5]]
        new_data[3] = "pk%i" % i
        print "\t".join(new_data)
