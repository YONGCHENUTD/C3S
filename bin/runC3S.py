#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 31 Aug 2017 11:35:46 AM
#
#         Module/Scripts Description
# 
# Copyright (c) 2016 Yunfei Wang <yfwang0405@gmail.com>
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  experimental
# @version: 1.0.0
# @author:  Yunfei Wang
# @contact: yfwang0405@gmail.com

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import pandas
import argparse


import c3s

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Gene Set Enrichment Analysis. Implemented according to Subramanian, A., et al. (2005). "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." Proc Natl Acad Sci U S A 102(43): 15545-15550.',add_help=False,epilog='dependency numpy, scipy, pandas')

    pr = p.add_argument_group('Required')
    pr.add_argument("-x","--genome",dest="genome",type=str,metavar="hg38", required=True, help="Bowtie2 built genome.")
    pr.add_argument("-1",dest="fq1",type=str,metavar='sample_R1.fastq.gz',required=True,help="Read 1 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.")
    pr.add_argument("-2",dest="fq2",type=str,metavar='sample_R2.fastq.gz',required=True,help="Read 2 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.")
    pr.add_argument("--prefix",dest="prefix",type=str,metavar='prefix',required=True,help="Prefix of result files.")

    po = p.add_argument_group('Optional')
    po.add_argument("--bait",dest="bait",type=str,metavar="chr11:5305934",default="chr11:5305934",help="Bait genomic locus. [Default=\"chr11:5305934\"]")
    po.add_argument("--extendsize",dest="extendsize",type=int,metavar="10000",default=100000,help="Length to be extended from bait regions. [Defaut=10000]")
    po.add_argument("--readlen",dest="readlen",type=int,metavar="36",default=36,help="Read length. [Default=36]")
    po.add_argument("--seed",dest="seed",type=int,metavar="1024",default=1024,help="Seed to generate random values. [Default=1024].")
    po.add_argument("--smooth-window",dest="smooth_window",type=str,metavar="100",default=100,help="Smooth window for peak size inference. [Default=100].")
#    po.add_argument("--max_gene",dest="maxgene",type=int,metavar="500",default=500,help="Maximum number of genes in gene set.")
    po.add_argument("-w",dest="wdir",type=str,metavar='"."',default=".",help="Working directory. [Default=\".\"].")
    po.add_argument("-p",dest='proc',type=int,metavar='10',default=10,help="Number of processes. [Default=10]")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args


# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    args = argParser()
    # Mapping reads to genome
    wdir = args.wdir+"/010ReadMapping"
    c3s.Utils.touchdir(wdir)
    wdir += "/"

    # 1st round of mapping
    c3s.Utils.touchtime("FIRST ROUND OF MAPPING ...")
    c3s.Utils.touchtime("MAPPING READ 1 ...")
    c3s.Tools.bowtie2_SE(args.genome,args.fq1,args.prefix+"_R1",proc=args.proc,wdir=wdir)
    c3s.Utils.touchtime("MAPPING READ 2 ...")
    c3s.Tools.bowtie2_SE(args.genome,args.fq2,args.prefix+"_R2",proc=args.proc,wdir=wdir)
    c3s.Utils.touchtime()

    # Split the reads by GATC sites and take the larger one
    c3s.Utils.touchtime("Split read by GATC sites ...")
    c3s.Algorithms.ParseGATCSites("{0}/{1}_R1_un.fastq.gz".format(wdir,args.prefix),"{0}/{1}_R1_split.fastq.gz".format(wdir,args.prefix))
    c3s.Algorithms.ParseGATCSites("{0}/{1}_R2_un.fastq.gz".format(wdir,args.prefix),"{0}/{1}_R2_split.fastq.gz".format(wdir,args.prefix))
    c3s.Utils.touchtime()

    # 2nd round of mapping
    c3s.Utils.touchtime("SECOND ROUND OF MAPPING ...")
    c3s.Utils.touchtime("MAPPING READ 1 ...")
    c3s.Tools.bowtie2_SE(args.genome,"{0}/{1}_R1_split.fastq.gz".format(wdir,args.prefix),args.prefix+"_R1_remap",min_qual=30,proc=args.proc,wdir=wdir)
    c3s.Utils.touchtime("MAPPING READ 2 ...")
    c3s.Tools.bowtie2_SE(args.genome,"{0}/{1}_R2_split.fastq.gz".format(wdir,args.prefix),args.prefix+"_R2_remap",min_qual=30,proc=args.proc,wdir=wdir)
    c3s.Utils.touchtime()

    # Fix mate pairs
    c3s.Utils.touchtime("Merge bam files and fix mate pairs ...")
    bams = [wdir+args.prefix+f for f in ["_R1.bam", "_R2.bam", "_R1_remap.bam", "_R2_remap.bam"]]
    tbffile = c3s.Algorithms.FixMatePairs(bams,wdir+args.prefix,args.proc)
    c3s.Utils.touchtime()

    # Infer peak characteristics from the bait region
    c3s.Utils.touchtime("Draw bait figures ...")
    peaksize = c3s.Plot.BaitStatsPlot(tbffile,args.bait,
                                      wdir+args.prefix+"_stats.pdf",
                                      extendsize=args.extendsize,
                                      readlen=args.readlen,
                                      smooth_window=args.smooth_window)
    c3s.Utils.touchtime()

    # Calculate intra- and inter-chrom interactions



