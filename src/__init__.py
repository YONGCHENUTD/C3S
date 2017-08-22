#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 22 Aug 2017 10:31:19 AM
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
import copy
from subprocess import call,PIPE

# non-built-in packages
import numpy
import pysam
import pandas

# ------------------------------------
# constants
# ------------------------------------

debug = True

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

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

class TabixFile(object):
    '''
    Class for Tabix file.
    '''
    def __init__(self,tbfile):
        self.closed = True
        self.infile = tbfile
        if not os.path.isfile(self.infile+".tbi"):
            self.infile = pysam.tabix_index(self.infile,seq_col=0,start_col=1,end_col=1,zerobased=True)
        self.fh = pysam.Tabixfile(self.infile)
        self.closed = False
    def __enter__(self):
        ''' Enter instance. '''
        return self
    def __exit__(self,etype,value,traceback):
        ''' Exit instance. '''
        self.close()
    def __del__(self):
        ''' On deletion. '''
        self.close()
    def close(self):
        ''' close file handle. '''
        if not self.closed:
            self.fh.close()
            self.closed = True
    def GetDepth(self,chrom,start,end,readlen=36):
        '''
        Get the depth of the corresponding region.
        Parameters:
            chrom: string
                chromosome
            start, end: int
                genomic coordinates.
            readlen: int
                read length
        Returns:
            depth: numpy.array
                read depth.
        '''
        depth = numpy.zeros(end-start)
        for item in self.fh.fetch(reference=chrom,start=start,end=end):
            tstart = int(item.split()[1]) - start
            tend = min(tstart+readlen,end-start)
            depth[tstart:tend] += 1
        return depth

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
    def ParseGATCSites(fqfile,outfile,overwrite=False):
        '''
        Parse GATC sites in unmapped reads. Split the reads into two parts, and keep the larger one.
        '''
        if not overwrite and os.path.isfile(outfile):
            Utils.touchtime("Output file exists: {0}, skipped ...".format(outfile))
            return
        total = 0
        cnt = 0
        with gzip.open(outfile,'wb') as ofh:
            for fq in IO.seqReader(fqfile,'fastq'):
                idx = fq.seq.find('GATC')  
                total += 1
                if idx == -1:
                    continue
                # report the larger one
                cnt += 1
                if idx< len(fq)/2-2: 
                    fq.seq = fq.seq[idx:]
                    fq.qual = fq.qual[idx:]            
                else:
                    fq.seq = fq.seq[:idx+4]
                    fq.qual = fq.qual[:idx+4]
                print >>ofh, fq
        Utils.touchtime("Read with GATC sites: {0} out {1} reads.".format(cnt,total))
    ParseGATCSites=staticmethod(ParseGATCSites)
    def _RmDup():
        '''
        Removed duplicate read pairs.
        '''
        pairs = {}
        curchr, curpos = "", 0
        reads = set()
        for line in sys.stdin:
            if line.startswith('@'):
                # get chromosomes
                if line[1:3].upper() == 'SQ':
                    pairs[line.split()[1].split(":")[1]] = {} # SN:chr1
                continue
            items = line.split('\t')
            items[3], items[7] = int(items[3]), int(items[7]) 
            if (items[6]== "=" and items[7]>items[3]) or items[2] < items[6]: # read1 < read2
                if items[2] == curchr and items[3] == curpos:
                    if items[6] == '=':
                        items[6] = curchr
                    reads.add("{0}:{1:0>10}".format(items[6],items[7])) # read pair with the same starts only count once
                else:
                    for read in reads:
                        chrom, pos = read.split(':')
                        pos = int(pos)
                        # read 1
                        pairs[curchr].setdefault(curpos,{})
                        pairs[curchr][curpos].setdefault(chrom,[])
                        pairs[curchr][curpos][chrom].append(pos)
                        # read 2
                        pairs[chrom].setdefault(pos,{})
                        pairs[chrom][pos].setdefault(curchr,[])
                        pairs[chrom][pos][curchr].append(curpos)
                    curchr, curpos = items[2], items[3]
                    reads = set()
        # last cycle
        for read in reads:
            chrom, pos = read.split(':')
            # read 1
            pairs[curchr].setdefault(curpos,{})
            pairs[curchr][curpos].setdefault(chrom,[])
            pairs[curchr][curpos][chrom].append(pos)
            # read 2
            pairs[chrom].setdefault(pos,{})
            pairs[chrom][pos].setdefault(curchr,[])
            pairs[chrom][pos][curchr].append(curpos)
        for curchr in sorted(pairs):
            for curpos in sorted(pairs[curchr]):
                for chrom in sorted(pairs[curchr][curpos]):
                    for pos in sorted(pairs[curchr][curpos][chrom]):
                        print "{0}\t{1}\t{2}\t{3}".format(curchr,curpos,chrom,pos)
    _RmDup=staticmethod(_RmDup)
    def FixMatePairs(bams,prefix,nproc=1,overwrite=False):
        '''
        Identify read pairs from bam files generated from two rounds of mapping.
        '''
        # check output file
        if not overwrite and os.path.isfile("{0}.pairs".format(prefix)):
            Utils.touchtime("Output file exists: {0}.pairs. Skipped ...".format(prefix))
            return
        # all on the fly
        # merge SE reads
        # fix mate pairs
        # select paired reads
        # sort by coordinates and print in SAM format
        # Remove duplicates, keep the left reads of the pairs
        cmd = '''samtools merge -nf -h {0} - {1} |samtools fixmate -pr - - |samtools view -h -f 1|samtools sort -@ {2} -O SAM -T {3} -|python -c "import c3s;c3s.Algorithms._RmDup()" >{3}.pairs'''.format(bams[0]," ".join(bams), max(nproc-3,1), prefix)
        Utils.touchtime(cmd)
        rc = call(cmd,shell=True)
        if rc != 0:
            raise ValueError("ERROR: FixMatePairs failed.")
        # creat tabix index
        Utils.touchtime("Build index for {0}.pairs".format(prefix))
        pairfile = pysam.tabix_index("{0}.pairs".format(prefix))
        Utils.touchtime("FixMatePairs finished ...")
        return pairfile
    FixMatePairs=staticmethod(FixMatePairs)
    def DetermineBinsize(depth,smooth_window=100):
        '''
        Determine the binsize of a region. The binsize is defined as the continous region around the bait (middle) position whose average window depth is larger than the average depth of the whole depth array.
        Paramegters:
            depth: numpy.array
                read depth
            smooth_window: int
                smooth window to determine the boundary.
        Returns:
            binsize: int
                size of the bait region.
        '''
        avgd = depth.mean()*smooth_window
        n, mid = len(depth), len(depth)/2
        left, right = 0, n
        # right boundary
        start = depth[mid]
        winsum = depth[mid:mid+smooth_window].sum()
        for i in range(mid+smooth_window,n):
            if winsum < avgd:
                right = i-1
                break
            winsum += depth[i] - start
            start = depth[i-smooth_window+1]
        # left boundary
        start = depth[mid-1]
        winsum = depth[mid-smooth_window:mid].sum()
        for i in range(mid-smooth_window-1,-1,-1):
            if winsum < avgd:
                left = i+1
                break
            winsum += depth[i] - start
            start = depth[i+smooth_window-1]
        return right-left
    DetermineBinsize=staticmethod(DetermineBinsize)

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
    def bowtie2_SE(genome, fqfile, prefix, proc=10, wdir=".", min_qual=0, overwrite=False):
        '''
        Bowtie version 2
        cmd: 
            >bowtie2 [options]* -x <bt2-idx> -U <r>} --un-gz <unmapped.fastq.gz |samtools view -Sb -q <min_qual> -F 4 - |samtools sort -n - -o <out.bam>
        The mapped reads are sorted by name, so that it would be easier when find corresponding pairs later.
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
        cmd = '''bowtie2 -x {genome} -U {fqfile} --un-gz {prefix}_un.fastq.gz {proc1} 2>{prefix}_bowtie2.log |\
samtools view -Sb {qual} -F 4 - |\
samtools sort -n {proc2} - -o {prefix}.bam 2>{prefix}_samtools.log
'''.format(genome=genome,
           fqfile=fqfile, 
           prefix=prefix, 
           proc1= "-p {0}".format(proc1) if proc1>1 else "",
           proc2="-@ {0}".format(proc2) if proc2>1 else "",
           qual= "-q {0}".format(min_qual) if min_qual else "")
        Utils.touchtime("Running command: {0}".format(cmd.replace("\\\n","")).rstrip())
        rc = call(cmd, shell=True)
        if rc != 0:
            os.remove(samfile)
        Utils.touchtime("Bowtie2 finishes ...")
    bowtie2_SE=staticmethod(bowtie2_SE)
 
# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example:"+sys.argv[0]+" genome fqfile prefix ")
    Tools.bowtie2_SE(sys.argv[1], sys.argv[2],sys.argv[3])

