#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 20 Jul 2017 07:05:45 AM
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

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class Algorithms(object):
    '''
    '''
    def test():
        pass


class Utils(object):
    '''
    Utilities.
    '''
    def cmdmustexist(cmd):
        ''' Test if commond exists. '''
        if call("type " + cmd, shell=True,  stdout=PIPE, stderr=PIPE) != 0:
            raise ValueError("ERROR: {0} cannot be found.".format(cmd))
    cmdmustexist=staticmethod(cmdmustexist)
    def mustexist(f):
        ''' check if a file exists, or raise an error. '''
        if not os.path.isfile(f):
            raise ValueError('ERROR: {0} does not exist.'.format(f))
    mustexist=staticmethod(mustexist)
    def touchtime(lstr,ofh=sys.stderr):
        ofh.write("# {0}: {1}\n".format(time.ctime(),lstr))
    touchtime=staticmethod(touchtime)

class Tools(object):
    '''
    '''
def bowtie2_SE(genome, fqfile, prefix, proc=10, overwrite=False):
        '''
        Bowtie version 2
        bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
        '''
        Utils.touchtime("Checking tools and input files.")
        Utils.cmdmustexist('bowtie2')
        Utils.cmdmustexist('samtools')
        Utils.touchtime("Starting bowtie2 ...")

        # check if samfile exists
        if os.path.isfile(samfile):
            if overwrite:
                os.remove(samfile)
            else:
                Utils.touchtime("Skipped: output file: {0} exists.\n".format(samfile))
                return

        # check if fastq file exists
        Utils.mustexist(fqfile)
        # check if genome exists
        Utils.mustexist(genome+".1.bt2")

        # run bowtie
        if proc >1:
            proc2 = int(proc/3.)
            proc1 = proc - proc2
        cmd = "bowtie2 -x {genome} -U {fqfile} -S --un {prefix}_un.fastq.gz {proc1} 2>{prefix}.log |samtools view -Sb -F 4 - |samtools sort {proc2}- {prefix} 2>{prefix}_samtools.log".format(genome=genome, fqfile=fqfile, prefix=prefix, proc1= "-p {0}".format(proc1) if proc1>1 else "",proc2="-@ {0}".format(proc2) if proc2>1 else "")
        Utils.touchtime("Running command: {0}\n".format(cmd))
        rc = call(cmd, shell=True)
        if rc != 0:
            os.remove(samfile)
        Utils.touchtime("Bowtie2 finishes ...")
    bowtie2=staticmethod(bowtie2)
 
# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example:"+sys.argv[0]+" genome fqfile prefix "
    bowtie2_SE(sys.argv[1], sys.argv[2],sys.argv[3])

