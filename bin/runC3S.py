#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 12 June 2020 10:48:53 PM
#
#         Module/Scripts Description
# 
# Copyright (c) 2020 Rowan University
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @version: 1.1.0
# @design: Yong Chen <chenyong@rowan.edu>
# @implementation: Yunfei Wang <yfwang0405@gmail.com>
# @corresponding author:  Yong Chen <chenyong@rowan.edu>

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
    p=argparse.ArgumentParser(description='C3S: model-based analysis and pipeline of dCas9 Capture-3C-Seq data.',add_help=False,epilog='dependency numpy, scipy, pandas, pysam, statsmodels')

    pr = p.add_argument_group('Required')
    pr.add_argument("-x","--genome",dest="genome",type=str,metavar="hg38", required=True, help="Bowtie2 built genome.")
    pr.add_argument("-1",dest="fq1",type=str,metavar='sample_R1.fastq.gz',nargs="+",required=True,help="Read 1 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.")
    pr.add_argument("-2",dest="fq2",type=str,metavar='sample_R2.fastq.gz',nargs="+",required=True,help="Read 2 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.")
    pr.add_argument("--prefix",dest="prefix",type=str,metavar='prefix',required=True,help="Prefix of result files.")
    pr.add_argument("--bait",dest="bait",type=str,metavar="chr11:5305934",required=True,help="Bait genomic locus. e.g. \"chr11:5305934\"")

    po = p.add_argument_group('Optional')
    po.add_argument("--extendsize",dest="extendsize",type=int,metavar="100000",default=100000,help="Length to be extended from bait regions. [Defaut=100000]")
    po.add_argument("--readlen",dest="readlen",type=int,metavar="36",default=36,help="Read length. [Default=36].")
    po.add_argument("--seed",dest="seed",type=int,metavar="1024",default=1024,help="Seed to generate random values. [Default=1024].")
    po.add_argument("--smooth-window",dest="smooth_window",type=int,metavar="100",default=100,help="Smooth window for peak size inference. [Default=100].")
    po.add_argument("--model-number",dest="nbins",type=int,metavar="11",default=11,help="Number of NB models used for testing. [Default=11].")
    po.add_argument("--nperm",dest="nperm",type=int,metavar="10000",default=10000,help="Number of permutatons. [Default=10000].")
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
    mappingdir = args.wdir+"/010ReadMapping"
    fq1, fq2 = ",".join(args.fq1), ",".join(args.fq2)
    mappingdir = c3s.Utils.touchdir(mappingdir)

    # 1st round of mapping
    c3s.Utils.touchtime("FIRST ROUND OF MAPPING ...")
    c3s.Utils.touchtime("MAPPING READ 1 ...")
    c3s.Algorithms.bowtie2_SE(args.genome,fq1,args.prefix+"_R1",proc=args.proc,wdir=mappingdir,min_qual=30)
    c3s.Utils.touchtime("MAPPING READ 2 ...")
    c3s.Algorithms.bowtie2_SE(args.genome,fq2,args.prefix+"_R2",proc=args.proc,wdir=mappingdir,min_qual=30)
    c3s.Utils.touchtime()

    # Split the reads by GATC sites and take the larger one
    c3s.Utils.touchtime("Split read by GATC sites ...")
    c3s.Algorithms.ParseGATCSites("{0}/{1}_R1_un.fastq.gz".format(mappingdir,args.prefix),"{0}/{1}_R1_split.fastq.gz".format(mappingdir,args.prefix))
    c3s.Algorithms.ParseGATCSites("{0}/{1}_R2_un.fastq.gz".format(mappingdir,args.prefix),"{0}/{1}_R2_split.fastq.gz".format(mappingdir,args.prefix))
    c3s.Utils.touchtime()

    # 2nd round of mapping
    c3s.Utils.touchtime("SECOND ROUND OF MAPPING ...")
    c3s.Utils.touchtime("MAPPING READ 1 ...")
    c3s.Algorithms.bowtie2_SE(args.genome,"{0}/{1}_R1_split.fastq.gz".format(mappingdir,args.prefix),args.prefix+"_R1_remap",min_qual=30,proc=args.proc,wdir=mappingdir)
    c3s.Utils.touchtime("MAPPING READ 2 ...")
    c3s.Algorithms.bowtie2_SE(args.genome,"{0}/{1}_R2_split.fastq.gz".format(mappingdir,args.prefix),args.prefix+"_R2_remap",min_qual=30,proc=args.proc,wdir=mappingdir)
    c3s.Utils.touchtime()

    # Fix mate pairs
    c3s.Utils.touchtime("Merge bam files and fix mate pairs ...")
    bams = [mappingdir+args.prefix+f for f in ["_R1.bam", "_R2.bam", "_R1_remap.bam", "_R2_remap.bam"]]
    tbffile = c3s.Algorithms.FixMatePairs(bams,mappingdir+args.prefix,args.proc)
    c3s.Utils.touchtime()

    # Infer peak characteristics from the bait region
    plotdir = args.wdir+"/020Plotting"
    plotdir = c3s.Utils.touchdir(plotdir)
    c3s.Utils.touchtime("Draw bait figures ...")
    tbf = c3s.TabixFile(tbffile)
    tbf.setChromSizes(bams[0])    
    tbf.BaitStatsPlot(args.bait,
                      plotdir+args.prefix+"_stats.pdf",
                      extendsize=args.extendsize,
                      readlen=args.readlen,
                      smooth_window=args.smooth_window)
    c3s.Utils.touchtime()

    # Calculate intra- and inter-chrom interactions
    modeldir = args.wdir+"/030Model"
    modeldir = c3s.Utils.touchdir(modeldir)
    c3s.Utils.touchtime("Permutation on intra-chromosomal interactions ...")
    ns, ps = tbf.GetIntraChromLinks(nperm=args.nperm,nbins=args.nbins)
    #for n,p in zip(ns,ps):
    #    print n,p
    c3s.Utils.touchtime("Permutation on inter-chromosomal interactions ...")
    n, p = tbf.GetInterChromLinks(nperm=args.nperm)
    
    #print n, p
    c3s.Utils.touchtime()

    # Calculate p values for intra- and inter-chrom interactions.
    c3s.Utils.touchtime("Calculate p values ...")
    tbf.InferBaitPval(modeldir+args.prefix)
    c3s.Utils.touchtime()

    # Ending
    c3s.Utils.touchtime("RunC3S finished successfully. Thank you for using C3S!")



