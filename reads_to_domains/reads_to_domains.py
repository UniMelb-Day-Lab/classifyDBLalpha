from mungo.fasta import FastaReader
from collections import defaultdict
from optparse import OptionParser



def main():
    parser = OptionParser()

    parser.add_option("", "--hmm", dest="hmmout",
        help="the nhmmOut file from running allocate_reads.py")

    parser.add_option("-o", "--out", dest="outfile",
        help="the location of an output file")

    (options, args) = parser.parse_args()

    top_hits = defaultdict(lambda: ("", 100.0))

    #Get top hits for each sequence
    for line in open(options.hmmout, 'rU'):
        if line[0]=="#": continue
        line=line.strip().split()
        read = line[0]
        E = float(line[4])
        if  E < top_hits[read][1]:
            top_hits[read] = (line[2], E)

    with open(options.outfile, 'w') as outfile:
     outfile.write("read,domain,e-value\n")
     for read in top_hits:
       outfile.write(read+","+top_hits[read][0]+","+str(top_hits[read][1])+"\n")

    return


if __name__ == '__main__':
    main()
