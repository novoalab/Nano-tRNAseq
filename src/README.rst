src
---

Installation
============

Nano-tRNAseq is written in Python3 and should work in most UNIX systems.

Make sure you install all programs listed below, before runnning the pipeline.
It's recommended to create virtual environment and use Python version 3.7.

.. code-block:: bash
		
   # create & activate new virtual environment
   mkdir -p ~/src/venv
   python3 -m venv ~/src/venv/Nano-tRNAseq
   source ~/src/venv/Nano-tRNAseq/bin/activate
   
   # install dependencies
   pip install numpy seaborn pandas mappy pysam

Note, this will install only Python dependencies.
We assume you have tools such as ``bwa mem`` and ``samtools``
already installed in your system.
If not, those can be easily installed with `conda <https://bioconda.github.io/>`_:

.. code-block:: bash
		
   conda install bwa samtools
   

Alignment
=========

Note, below we'll use only 10% of randomly selected reads from replicate 1.
This will align tRNA reads from ``guppy/*/``:

.. code-block:: bash

   cd ~/src/Nano-tRNAseq/test
   
   ref=~/src/Nano-tRNAseq/ref/yeast.tRNA.ext.fa
   bwa index $ref

   params="-W 13 -k 6 -x ont2d -T20"
   mkdir -p bwamem
   for f in guppy/*/; do
     s=`basename $f`
     echo `date` $f $s
     bwa mem -t6 $params $ref <(zcat $f/*.fastq.gz|sed 's/U/T/g') 2> /dev/null | samtools sort --write-index -o bwamem/$s.bam
   done; date

Note, you'll need to basecall the reads before.


Estimate expression
===================

The expression of tRNAs can be estimated as follows:

.. code-block:: bash
		
   ~/src/Nano-tRNAseq/src/get_counts.py -o rep1.counts.tsv \
     -i bwamem/{wt,pus4ko,h2o2}.bam > rep1.stats.tsv

      
This will produce an output file ``rep1.counts.tsv`` similar to this:

.. code-block:: bash

    ref     wt      pus4ko  h2o2
    Ala-AGC 1175    1103    847
    Ala-TGC 472     314     419
    Arg-ACG 1547    1342    1019
    ...
    iMet-CAT        81      182     154
    RDN5    8231    4936    3051
    RDN58   1400    239     532
    oligo3  848     474     582
    oligo5  36      27      20
    antisense       45      27      11
    not_unique      16131   8050    7118
    *       44947   32046   27221


Adittionally, it'll report a statistics table to ``rep1.stats.tsv`` like this one

.. code-block:: bash

    sample  reads   aligned %       unique tRNA     %       antisense       %       oligo3  %       oligo5  %
    wt      116139  71192   61.3    44501   62.51   45      0.1     848     1.19    36      0.05
    pus4ko  78583   46537   59.22   32784   70.45   27      0.08    474     1.02    27      0.06
    h2o2    65275   38054   58.3    26740   70.27   11      0.04    582     1.53    20      0.05


Finally, you can plot scatterplot of expression for all
(or for selected ``--sample``) using: 

.. code-block:: bash

    ~/src/Nano-tRNAseq/src/plot_scatter.py -i rep1.counts.tsv -o rep1.scatter.pdf

	   
Modification detection
======================

First, you'll need to generate difference in sum of basecalling errors
between WT and some other sample(s) as follows:

.. code-block:: bash
		
   ~/src/Nano-tRNAseq/src/get_sum_err.py -o heatmap/rep1 -f $ref \
     -i bwamem/{wt,pus4ko,h2o2}.bam

Optionally, you can add modification annotation using
``-b ~/src/Nano-tRNAseq/ref/yeast.tRNA.modomics.tsv.bed``. 
     
Then, the heatmaps can be plotted using:

.. code-block:: bash
		
   ~/src/Nano-tRNAseq/src/plot_heatmap.py -i heatmap/rep1.err_diff.tsv.gz \
     -f $ref -a ~/src/Nano-tRNAseq/ref/yeast.tRNA.modomics.aln.fa

