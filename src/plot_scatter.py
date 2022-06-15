#!/usr/bin/env python3
desc="""Plot scatterplot
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 27/4/2022
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

def plot_scatter(outfn, tsv, samples=[]):
    """Plot scatterplot"""
    print("Loading data...")
    if len(samples)==1:
        sys.stderr.write("You have to specify 0 or 2+ samples\n")
        sys.exit(1)
    # load counts
    df = pd.read_csv(tsv, sep="\t")
    # and get only tRNAs
    sel = np.array([True if "-" in ref and ref!="not-unique" else False
                    for ref in df["ref"]])
    df2 = df[sel].copy()
    # add Aminoacid
    df2["Aminoacid"] = df2['ref'].str.slice(-7, -4)
    #
    title = "Spearman's Rho=%.3f"
    cols = df2.columns[1:-1]
    if samples:
        for s in samples:
            if s not in cols:
                sys.stderr.write("[ERROR] %s not in %s\n"%(s, tsv))
                sys.exit(1)
        cols = samples
    n_cols = len(cols)-1
    print("Plotting %s scatterplot(s) for: %s ..."%(n_cols, ", ".join(cols)))
    fig, axes = plt.subplots(1, n_cols, figsize=(7*n_cols, 7))
    if n_cols==1: axes = [axes, ]
    x = cols[0]
    legend = True
    for ax, y in zip(axes, cols[1:]):
        sys.stderr.write(" %s         \r"%y)
        sns.scatterplot(data=df2, x=x, y=y, ax=ax, legend=legend, hue="Aminoacid")
        ax.set_xscale('log'); ax.set_yscale('log')
        rho = df2[[x, y]].corr(method='spearman').to_numpy()[1,0]
        ax.set_title(title%rho)
        if legend:
            ax.legend(ncol=2)
        # plot legend only on the first panel
        legend = False
    sns.despine()
    fig.savefig(outfn)
    print("Figure saved as %s ."%outfn)

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b') 
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--input", default="counts.tsv", help="input file")
    parser.add_argument("-o", "--out", default="scatter.pdf",
                        help="output name [%(default)s]")
    parser.add_argument("-s", "--sample", nargs="+", default=[], 
                        help="sample names [use all]")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    plot_scatter(o.out, o.input, o.sample)
    
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
