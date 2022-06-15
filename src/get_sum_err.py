#!/usr/bin/env python3
desc="""Generate basecalling error differeces.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 25/4/2022
"""

import os
import sys
import pysam
import numpy as np
import pandas as pd
from datetime import datetime

### modified from https://github.com/lpryszcz/REDiscover
alphabet = "ACGTdi" #se
base2name = {"A": "A", "C": "C", "G": "G", "T": "T",
             "i": "insertions", "d": "deletions"}
base2index = {b: i for i, b in enumerate(alphabet)}
for i, b in enumerate(alphabet.lower()): base2index[b] = i
# code N as A
base2index["N"] = 0

# CIGAR operations
"""Op BAM Description +1Q +1R
M 0 alignment match (can be a sequence match or mismatch) yes yes
I 1 insertion to the reference yes no
D 2 deletion from the reference no yes
N 3 skipped region from the reference no yes
S 4 soft clipping (clipped sequences present in SEQ) yes no
H 5 hard clipping (clipped sequences NOT present in SEQ) no no
P 6 padding (silent deletion from padded reference) no no
= 7 sequence match yes yes
X 8 sequence mismatch yes yes
"""
def _match(refi, readi, bases): return refi+bases, readi+bases, True
def _insertion(refi, readi, bases): return refi, readi+bases, False
def _deletion(refi, readi, bases): return refi+bases, readi, False
def _skip(refi, readi, bases): return refi, readi, False
code2function = {0: _match, 7: _match, 8: _match, 1: _insertion, 6: _insertion,
                 2: _deletion, 3: _deletion, 4: _insertion, 5: _skip}

def store_blocks(a, start, end, baseq, i, calls):
    """Store base calls from aligned blocks. INDEL aware."""
    # process read blocks and store bases counts and indels as well
    readi, refi = 0, a.pos
    for ci, (code, bases) in enumerate(a.cigar):
        prefi, preadi = refi, readi
        refi, readi, data = code2function[code](refi, readi, bases)
        # skip if current before start
        if refi<=start: continue
        # typical alignment part
        if data:
            if prefi<start:
                bases -= start-prefi
                preadi += start-prefi
                prefi = start
            if refi>end: bases -= refi-end
            if bases<1: break
            for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases])):
                if q>=baseq and b in base2index:
                    calls[prefi-start+ii, i, base2index[b]] += 1
        elif start<prefi<end:
            # insertion
            if code==1: calls[prefi-start, i, 5] += 1
            # deletion
            elif code==2: calls[prefi-start, i, 4] += 1
    return calls

def bam2calls(sam, ref, start=0, end=None, mapq=1, baseq=0):
    """Return 3D array of basecalls from BAM file, as follows:
    - 1D positions from start to end of the ref
    - 2D sense and antisense strand
    - 3D base counts for ACGTidse
    """
    if not end: 
        ref2len = {r: l for r, l in zip(sam.references, sam.lengths)}
        end = ref2len[ref]
    # position, strand, ACGTid
    calls = np.zeros((end-start+1, 2, len(alphabet)), dtype="int64")
    # stop if ref not in sam file
    for a in sam.fetch(ref, start, end):
        if a.mapq<mapq or a.is_secondary or a.is_qcfail: continue
        # get transcript strand
        i = 0 # for +/for i == 0; for -/rev i==1
        if a.is_reverse: i = 1
        # store alignment blocks
        calls = store_blocks(a, start, end, baseq, i, calls)
    return calls
###

def mapped_uniquely(a, mapq=1):
    if a.is_secondary or a.is_supplementary or a.mapq<mapq:
        return
    return True

def get_errors(outfn, mods, faidx, refs, bams, samples, adapter5len, ignore_indels=False):
    """Return df with errors"""
    err = [(ref, p) for ref in refs for p in range(1, len(faidx.fetch(ref))+1) #-adapter5len
           if "-" in ref]; err
    err = pd.DataFrame(data=err, columns=["chrom", "pos"]); err
    # add mods
    new_cols = ["name", ] # "Canonical position"
    for c in new_cols: err[c] = ""
    for ref in refs:
        pos = mods.loc[mods.chrom==ref, "end"]+adapter5len
        sel = (err.chrom==ref)&(err.pos.isin(pos))
        for c in new_cols:
            err.loc[sel, c] = mods.loc[mods.chrom==ref, c].to_list()
    
    cols = ["%s %s"%(s, base2name[b]) for s in samples for b in alphabet]; cols
    cols0 = ["%s %s"%(s, e) for s in samples for e in ["sum_err", ]]; cols0
    #refs = mods.chrom.unique()
    sams = [pysam.AlignmentFile(bam) for bam in bams]
    for ri, ref in enumerate(refs, 1):
        if ref not in sams[0].references: continue
        sys.stderr.write("%s / %s %s\r"%(ri, len(refs), ref))
        freq = []
        for sam in sams:
            calls = bam2calls(sam, ref)
            # get call freq - we use only ACGT and del for coverage, but not insertions
            f = (calls[:, 0, :].T / calls[:, 0, :5].sum(axis=1)).T
            # store only sense counts
            freq.append(f)
        # get array
        freq = np.array(freq)
        # get positions with mods
        pos = err.loc[err.chrom==ref].pos.to_numpy()-1#+adapter5len
        # calculate and store sum of errors
        seq = faidx.fetch(ref)
        idx = [base2index[seq[p]] for p in pos]; idx
        if ignore_indels:
            sum_err = [(1-freq[:, pos[i], idx[i]]) / freq[:, pos[i], :4].sum(axis=1) 
                       for i in range(len(pos))]
        else:
            sum_err = [(1-freq[:, pos[i], idx[i]]) / freq[:, pos[i]].sum(axis=1) 
                       for i in range(len(pos))]
        # store rounded .123 sum_err
        err.loc[err.chrom==ref, cols0] = np.round(sum_err, 3)
        # stack samples and round to .123
        freq_stacked = np.round(np.hstack(freq), 3)
        # process references with annotated mods - estimate mods freq
        err.loc[err.chrom==ref, cols] = freq_stacked[pos]
    # save tsv
    err.dropna(thresh=len(cols)).to_csv(outfn, sep="\t", index=False)
    return err

def compare_to_wt(err, samples):
    """Return df with difference between first sample the the other samples"""
    cols0 = list(err.columns[:3])
    se_cols = ['%s sum_err'%s for s in samples]
    # compare everything to wt_RT
    err["wt"] = err[se_cols[0]]
    cols = ["wt", ]
    for c in se_cols[1:]: 
        s = c.split()[0]
        cols.append(s)
        err[s] = err[c]-err[cols[0]]
    return err[cols0+cols]

def get_sum_err(bams, samples, bed, fasta, fname, adapter5len=24, adapter3len=30, 
                sort_idx=0, ignore_indels=False, overwrite=0):
    """Load alignments, modifictions and errors from BAM and store sum_err heatmap"""
    base_freq_fn = fname+".base_freq.tsv.gz"
    err_diff_fn = fname+".err_diff.tsv.gz"
    if os.path.isfile(err_diff_fn) and not overwrite:
        print("File exists: %s"%err_diff_fn)
        return err_diff_fn
    # get samples from file names
    if not samples:
        samples = [os.path.basename(fn)[:-4] for fn in bams]
    # get FastA
    faidx = pysam.FastaFile(fasta)
    refs = [r for r in faidx.references if "-" in r]
    # get refs sorted by decreasing number of uniquely mapped reads - use first bam
    sam = pysam.AlignmentFile(bams[sort_idx])
    refs_sorted = sorted(refs, reverse=True, key=lambda r: 
                         sam.count(r, read_callback=mapped_uniquely))    
    # create outdir
    if os.path.sep in fname:
        outdir = os.path.dirname(fname)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    # load mods
    cols = ["chrom", "start", "end", "name", "score", "strand", "s", "e", "rgb"]
    if bed:
        print("Loading modifications...")
        df = pd.read_csv(bed, sep="\t", names=cols)
        # mod names
        names = {n for n in df.name.unique() if not n.startswith("oligo")}; 
        mods = df.loc[df.name.isin(names), cols[:4]].copy()
        # -24 for start and end!!!
        mods.start -= adapter5len
        mods.end -= adapter5len
    else:
        mods = pd.DataFrame(data=[], columns=cols)
    # get errors
    print("Loading basecalling errors...")
    err = get_errors(base_freq_fn, mods, faidx, refs_sorted, bams, samples,
                     adapter5len, ignore_indels)
    print("Saving sum of error differences to %s ..."%err_diff_fn)
    err_diff = compare_to_wt(err, samples)
    # and write tsv file
    err_diff.to_csv(err_diff_fn, sep="\t", index=False)
    return err_diff_fn
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b') 
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")
    parser.add_argument("-i", "--bam", nargs="+", help="input BAM files")
    parser.add_argument("-s", "--sample", nargs="+", help="sample names")
    parser.add_argument("-o", "--out", default="heatmap/rep0",
                        help="output name [%(default)s]")
    parser.add_argument("-b", "--bed", default="", help="tRNA mods in BED format")
    parser.add_argument("-f", "--fasta", required=True, help="tRNA FastA file")
    parser.add_argument("--len5", default=24, type=int, help="oligo5 length [%(default)s]")
    parser.add_argument("--len3", default=30, type=int, help="oligo3 length [%(default)s]")
    parser.add_argument("--sort_idx", default=0, type=int, 
                        help="sample used for reference sorting [%(default)s]")
    parser.add_argument("--ignore_indels", action="store_true",
                        help="ignore indels for sum of errors calculation")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    get_sum_err(o.bam, o.sample, o.bed, o.fasta, o.out,
                o.len5, o.len3, o.sort_idx, o.ignore_indels)
    
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
