#!/usr/bin/env python3
desc="""Report read counts for references

TODO:
- report reads that aligned not-uniquely on single isoacceptor-anticodon
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 27/4/2022
"""

import os
import sys
import pysam
import pandas as pd
from datetime import datetime

def get_ref2count(bam, mapq=1):
    """Return per-reference read counts"""
    sam = pysam.AlignmentFile(bam)
    ref2c = {r: 0 for r in [*sam.references, "antisense", "not_unique", "*"]}
    for a in sam:
        # skip secondary and supplementary alignments
        if a.is_secondary or a.is_supplementary: continue
        # minimap2 reports None instead of *
        if a.reference_name: 
            if a.mapq>=mapq: 
                # count antisense alignments - those are false positives!
                if a.is_reverse: ref = "antisense"
                else: ref = a.reference_name
            # count not-uniquely aligned reads on oligo
            elif a.reference_name.startswith("oligo"):
                ref = a.reference_name
            else:
                ref = "not_unique"
        else:
            ref = "*"
        # store
        ref2c[ref] += 1
    return ref2c

def get_counts(outfn, bams, samples, mapq=1, log=sys.stderr):
    """Report per-reference read counts"""
    if not samples:
        # path/to/sample.bam > sample
        samples = [os.path.basename(bam[:-4]) for bam in bams]
    elif len(bams)!=len(samples):
        sys.stderr.write("no. of bams and samples have to match!\n")
        sys.exit(1)
    # load counts
    log.write("Loading read counts...\n")
    refs = []
    line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
    print(line%("sample", "reads", "aligned", "%", "unique tRNA", "%", "antisense", "%", 
                "oligo3", "%", "oligo5", "%"))
    percent = lambda x, y: round(100*x/y, 2) if y else 0
    for i, (bam, sample) in enumerate(zip(bams, samples)):
        log.write(" %s / %s %s   \r"%(i+1, len(bams), sample))
        ref2c = get_ref2count(bam, mapq)
        if not refs:
            refs = [ref for ref in ref2c]
            df = pd.DataFrame(data={"chrom": refs})
        df[sample] = [ref2c[ref] for ref in refs]
        # report stats
        reads = sum(ref2c.values())
        # here we assume that only tRNAs have `-` ie Ala-AGC
        unique = sum([v for r, v in ref2c.items() if "-" in r])
        antisense = ref2c["antisense"]
        unaligned = ref2c["*"]
        aligned = reads - unaligned
        oligo3, oligo5 = [ref2c[k] if k in ref2c else 0 for k in ("oligo3", "oligo5")]
        print(line%(sample, reads, aligned, percent(aligned, reads), 
                    unique, percent(unique, aligned), antisense, percent(antisense, unique), 
                    oligo3, percent(oligo3, aligned), oligo5, percent(oligo5, aligned)))        
    # sum counts for given isoacceptor-anticodon
    log.write("Summing counts for each isoacceptor-anticodon pair...\n")
    df["ref"] = ["-".join(r.split("-")[:2]) for r in df["chrom"]]
    df = df[df.columns[1:]].groupby("ref", sort=False).sum()
    # store counts
    log.write("Saving counts to %s ...\n"%outfn)
    df.to_csv(outfn, sep="\t")
    return df

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b') 
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--bam", nargs="+", help="input BAM files")
    parser.add_argument("-o", "--out", default="counts.tsv",
                        help="output name [%(default)s]")
    parser.add_argument("-s", "--sample", nargs="+", 
                        help="sample names [use file names]")
    parser.add_argument("-m", "--mapq", default=1, type=int, 
                        help="min. mapping quality [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    get_counts(o.out, o.bam, o.sample, o.mapq)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
