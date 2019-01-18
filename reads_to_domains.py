from mungo.fasta import FastaReader
from collections import defaultdict

top_hits = defaultdict(lambda: ("", 100.0))

#Get top hits for each sequence
for line in open("InsertFilename_nhmmOut.txt", 'rU'):
 if line[0]=="#": continue
 line=line.strip().split()
 read = line[0]
 E = float(line[4])
 if  E < top_hits[read][1]:
   top_hits[read] = (line[2], E)

with open("reads_to_domains.csv", 'w') as outfile:
 outfile.write("read,domain,e-value\n")
 for read in top_hits:
   outfile.write(read+","+top_hits[read][0]+","+str(top_hits[read][1])+"\n")
