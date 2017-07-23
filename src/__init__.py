#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 23 Jul 2017 02:10:53 AM
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
import gzip
import time
import errno
import pandas

from subprocess import call,PIPE

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

import copy
class Fasta(object):
    '''Fasta format.'''
    def __init__(self,name, seq, description=''):
        '''Initiate the fasta record.'''
        self.id=name
        self.seq = seq
        self.description=description
    def __len__(self):
        '''get the length of sequence.'''
        return len(self.seq)
    def length(self):
        '''get the length of sequence.'''
        return len(self.seq)
    def __str__(self):
        '''String for output of Fasta.'''
        return ">{0}{1}\n{2}".format(self.id,(self.description and ' '+self.description or ''),self.seq)
    
class Fastq(Fasta):
    ''' Fastq format. '''
    def __init__(self, name, seq, qual = '', description = ''):
        ''' Initiation. '''
        Fasta.__init__(self, name, seq, description)
        self.qual = qual
    def __str__(self):
        ''' String for output of Fastq. '''
        return "@{0}\n{1}\n+\n{2}".format(self.id,self.seq, len(self.qual)==self.length() and self.qual or ''.join(['I' for i in xrange(self.length())]))
    
class IO(object):
    def mopen(infile,mode='r'):
        ''' Open file with common or gzip types.'''
        if infile.endswith(".gz"):
            if 'b' not in mode:
                mode += 'b'
            return gzip.open(infile, mode)
        else:
            return open(infile,mode)
    mopen=staticmethod(mopen)
    def seqReader(infile,ftype='fasta'):
        '''Read sequence files.'''
        ftype=ftype.lower()
        # Read lines
        with IO.mopen(infile) as fh:
            if ftype=='fasta':
                line = fh.next()
                if line[0] != ">":
                    raise ValueError("Records in Fasta files should start with '>' character")
                line = line.lstrip('>').rstrip().replace('\t',' ').split(' ')
                name = line[0]
                desc = ' '.join(line[1:])
                seq = ''
                while True:
                    try:
                        line = fh.next()
                    except:
                        if seq != '':
                            yield Fasta(name, seq, desc)
                        raise StopIteration
                    if line[0] != ">":
                        seq += line.rstrip()
                    else:
                        yield Fasta(name, seq, desc)
                        line = line.lstrip('>').rstrip().replace('\t',' ').split(' ')
                        name = line[0]
                        desc = ' '.join(line[1:])
                        seq = ''
            elif ftype=='fastq':
                while True:
                    try:
                        fid=fh.next().rstrip().lstrip('@')
                        seq=fh.next().rstrip()
                        fh.next()
                        qual = fh.next().rstrip()
                        yield Fastq(fid,seq,qual)
                    except:
                        raise StopIteration
            else:
                raise TypeError(ftype+" format is not supported.")
            assert False, "Do not reach this line."
    seqReader=staticmethod(seqReader)

class Algorithms(object):
    '''
    '''
    def ParseGATCSites(fqfile,outfile):
        '''
        '''
        with gzip.open(outfile,'wb') as ofh:
            for fq in IO.seqReader(fqfile,'fastq'):
                idx = fq.seq.find('GATC')  
                if idx == -1:
                    continue
                # report the larger one
                if idx< len(fq)/2-2: 
                    fq.seq = fq.seq[idx:]
                    fq.qual = fq.qual[idx:]            
                else:
                    fq.seq = fq.seq[:idx+4]
                    fq.qual = fq.qual[:idx+4]
                print >>ofh, fq
    ParseGATCSites=staticmethod(ParseGATCSites)

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
    def touchtime(lstr="",ofh=sys.stderr):
        if lstr == "": # empty line as seperator
            ofh.write("\n")
        else:
            ofh.write("# {0}: {1}\n".format(time.ctime(),lstr))
    touchtime=staticmethod(touchtime)
    def touchdir(wdir="."):
        if not os.path.exists(wdir):
            try:
                os.makedirs(wdir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
    touchdir=staticmethod(touchdir)

class Tools(object):
    '''
    '''
    def bowtie2_SE(genome, fqfile, prefix, proc=10, wdir=".", overwrite=False):
        '''
        Bowtie version 2
        bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
        '''
        Utils.touchtime("Checking tools and input files.")
        Utils.cmdmustexist('bowtie2')
        Utils.cmdmustexist('samtools')
        Utils.touchtime("Starting bowtie2 ...")
        
        Utils.touchdir(wdir)
        wdir += "/"
        prefix = wdir+prefix
        # check if samfile exists
        samfile = prefix+".bam"
        if os.path.isfile(samfile):
            if overwrite:
                os.remove(samfile)
            else:
                Utils.touchtime("Skipped: output file: {0} exists.".format(samfile))
                return

        # check if fastq file exists
        Utils.mustexist(fqfile)
        # check if genome exists
        Utils.mustexist(genome+".1.bt2")

        # run bowtie
        if proc >1:
            proc2 = int(proc/3.)
            proc1 = proc - proc2
        cmd = "bowtie2 -x {genome} -U {fqfile} --un-gz {prefix}_un.fastq.gz {proc1} 2>{prefix}_bowtie2.log |samtools view -Sb -F 4 - |samtools sort {proc2} - -o {prefix}.bam 2>{prefix}_samtools.log".format(genome=genome, fqfile=fqfile, prefix=prefix, proc1= "-p {0}".format(proc1) if proc1>1 else "",proc2="-@ {0}".format(proc2) if proc2>1 else "")
        Utils.touchtime("Running command: {0}".format(cmd))
        rc = call(cmd, shell=True)
        if rc != 0:
            os.remove(samfile)
        else:
            Utils.touchtime("Build BAM index.")
            rc = call("samtools index {0}.bam".format(prefix),shell=True)
        Utils.touchtime("Bowtie2 finishes ...")
    bowtie2_SE=staticmethod(bowtie2_SE)
 
# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example:"+sys.argv[0]+" genome fqfile prefix ")
    Tools.bowtie2_SE(sys.argv[1], sys.argv[2],sys.argv[3])

