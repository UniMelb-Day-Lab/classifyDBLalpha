# Allocate NGS reads to protein domains
This python program wraps Uproc and HMMER3 allowing for the fast identification of NGS reads that belong to HMMER PHMMs. It was designed to be used with the protein domains of the Plasmodium falciparum VAR gene family but could be used for other domain models. By leverageing the speed of Uproc and the sensitivity of HMMER we are able to quickly annotate longer reads to smaller domain profiles. The code has been moved over from https://github.com/PapenfussLab/reads_to_domains.


### allocate_reads.py
```
Options:
  -h, --help            show this help message and exit
  -r READ1, --read1=READ1
                        the first fastq file
  -R READ2, --read2=READ2
                        the second fastq file
  -E EVALUE, --evalue=EVALUE
                        the E value to be fed to HMMER
  -o OUTDIR, --outdir=OUTDIR
                        the output directory for files
```

### reads_to_domains.py
```
Options:
  -h, --help            show this help message and exit
  --hmm=HMMOUT          the nhmmOut file from running allocate_reads.py
  -o OUTFILE, --out=OUTFILE
                        the location of an output file
```
