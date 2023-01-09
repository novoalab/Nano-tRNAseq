#!/usr/bin/env python3
desc="""Plot heatmaps of basecalling error differeces.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 25/4/2022
"""

import os
import sys
import pysam
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from datetime import datetime

def get_ticks(alnidx, ref2idx, mods, adapter5len):
    """Get ticks for the plot"""
    # for old yeast algs add 0,!
    ticks = "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,17a,18,19,20,20a,20b,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,V1,V2,V3,V4,V5,V11,V12,V13,V14,V15,V16,V17,V21,V22,V23,V24,V25,V26,V27,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76"
    ticks +=  ",C,C,A"
    shift = -1 if ticks.startswith("0") else 0
    ticks = ticks.split(',')
    aligned_mods = [[] for i in range(len(ticks))]
    if "name" in mods.columns:
        for ref, idx in ref2idx.items():
            pos2idx = {p+adapter5len: i for p, i in enumerate(idx)}
            for i, r in mods[mods.chrom==ref].iterrows():
                pos, name = r["pos"], r["name"]
                if pd.isna(name) or name.startswith("oligo"):
                    continue
                aligned_mods[pos2idx[pos]+shift].append(name)
    # combine ticks with mods
    ticks = ["%s %s"%(" ".join(set(m)), p) if m else p
             for p, m in zip(ticks, aligned_mods)]
    ticks = np.array(ticks)
    return ticks

def save_heatmaps(ticks, err_diff, ref2idx, fname, ext, ref2len, 
                  maxv=0.8, adapter5len=24):
    """Generate and save heatmap"""
    # set color pallete
    cmap = sns.color_palette("vlag", as_cmap=True) # vlag coolwarm    
    # get only refs that are aligned
    #_refs = [r for r in refs_sorted if r in ref2idx]; len(_refs)
    # err_diff from sum of errors
    if err_diff.columns[0] == "reference":
        df = err_diff.copy().fillna(0)
        df["WT"] = df[df.columns[-2:]].diff(axis=1)[df.columns[-1]]
        df = df.rename(columns={"reference": "chrom", "position": "pos"})
        err_diff = df[["chrom", "pos", "strand", "KO freq", "WT"]].copy()
        title = "mod_freq change %s vs KO"
    else:
        title = "sum_err change {} vs %s".format(err_diff.columns[3])
    _refs = [r for r in err_diff["chrom"].unique() if r in ref2idx]; len(_refs)
    coi = err_diff.columns[3:]; coi
    for c in coi[1:]:
        d = np.empty((len(_refs), len(ticks)))
        d[:] = np.nan
        for i, r in enumerate(_refs):
            idx = np.arange(ref2len[r], dtype=int)
            _data = err_diff.loc[(err_diff.chrom==r)&(err_diff.pos.isin(idx+adapter5len+1)), c]
            if len(ref2idx[r])!=len(_data):
                mssg = "Skipped %s due to length difference: %s!=%s"
                print(mssg%(r, len(ref2idx[r]), len(_data)))
                continue
            d[i, ref2idx[r]] = _data
        # get df
        amods = pd.DataFrame(data=d, columns=ticks)
        amods["chrom"] = _refs
        amods = amods.set_index("chrom")
        # drop rows & columns with only nans
        for axis in (0, 1):
            amods.dropna(axis=axis, how="all", inplace=True)
        # plot
        fig, ax = plt.subplots(figsize=(25, len(_refs)/3)) # for human we use /4
        g = sns.heatmap(amods, ax=ax, cmap=cmap, vmin=-maxv, vmax=maxv)
        ax.tick_params(axis="x", labelrotation=90)
        ax.set_title(title%c)
        # plot missing values in grey
        g.set_facecolor(".85") # 'xkcd:grey'
        fig.tight_layout()
        fig.savefig("%s.%s.%s"%(fname, c, ext))

def plot_heatmap(fname, fasta, alnfn, maxv, adapter5len=24, adapter3len=30, ext="pdf"):
    """Load alignments, modifictions and errors from BAM and store sum_err heatmap"""
    # load sum of errors
    err_diff = pd.read_csv(fname, sep="\t")
    # create outdir
    if os.path.sep in fname:
        outdir = os.path.dirname(fname)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    # get aligned positions
    print("Loading alignments...")
    faidx = pysam.FastaFile(fasta)
    alnidx = pysam.FastaFile(alnfn)
    ref2idx = {r: np.array([i for i, b in enumerate(alnidx[r]) if b!="-"])
               for r in alnidx.references if r in faidx}
    # get ticks
    ticks = get_ticks(alnidx, ref2idx, err_diff, adapter5len)
    # get ref lenths
    ref2len = {r: l - adapter3len - adapter5len
               for r, l in zip(faidx.references, faidx.lengths)}
    # generate & save heatmaps
    print("Plotting...")
    save_heatmaps(ticks, err_diff, ref2idx, fname, ext, ref2len, maxv, adapter5len)
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b') 
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")
    parser.add_argument("-i", "--tsv", required=True, help="input tsv file")
    parser.add_argument("-e", "--ext", default="pdf", choices=("pdf", "png", "jpg"), 
                        help="output format [%(default)s]")
    parser.add_argument("-f", "--fasta", required=True, help="tRNA FastA file")
    parser.add_argument("-a", "--aln", required=True, help="tRNA aligments")
    parser.add_argument("--maxv", default=0.8, type=float, 
                        help="max value in the heatmap [%(default)s]")
    parser.add_argument("--len5", default=24, type=int, help="oligo5 length [%(default)s]")
    parser.add_argument("--len3", default=30, type=int, help="oligo3 length [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    plot_heatmap(o.tsv, o.fasta, o.aln, o.maxv, o.len5, o.len3, o.ext)
    
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
