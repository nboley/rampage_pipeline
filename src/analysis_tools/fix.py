import sys
with open(sys.argv[1]) as fp:
    for line in fp:
        print "chr"+line,
