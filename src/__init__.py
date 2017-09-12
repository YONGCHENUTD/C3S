#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 12 Sep 2017 10:04:55 AM
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
import matplotlib
from scipy.stats import nbinom
from statsmodels.base.model import GenericLikelihoodModel

# ------------------------------------
# constants
# ------------------------------------

debug = True
matplotlib.use('Agg')

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
    def setChromSizes(self,bamfile):
        '''
        Set chromosome names and sizes.
        '''
        with pysam.Samfile(bamfile) as sam:
            self.chroms, self.sizes = sam.references,sam.lengths
    def BaitStatsPlot(self,bait,outfile,extendsize=100000,readlen=36,smooth_window=100):
        '''
        '''
        # Initiation
        Utils.touchtime("Initialize tabix file ...")
        bait_chrom,bait_pos = bait.split(':')
        bait_pos = int(bait_pos)
        self.bait_chrom, self.bait_pos = bait_chrom, bait_pos
        Utils.touchtime("Calculate peak region ...")
        start, end = bait_pos-extendsize, bait_pos+extendsize
        depth = numpy.zeros(end-start)
        for item in self.fh.fetch(reference=bait_chrom,start=start,end=end):
            items = item.split()
            pos, ochrom, opos = int(items[1]), items[2], int(items[3])
            tstart = int(item.split()[1]) - start
            tend = min(tstart+readlen,end-start)
            depth[tstart:tend] += 1

        # calculate the peaksize
        left, right = Algorithms.DeterminePeakSize(depth,smooth_window)
        sdepth = depth[left:right]
        Utils.touchtime("Peak size inferred from the bait region: {0}".format(right-left))
        left, right = left+start, right+start
        self.left, self.right = left, right
        self.peaksize = right-left

        # fetch bait peak links
        Utils.touchtime("Count bait related links ...")
        targets = {chrom:[] for chrom in self.chroms}
        for item in self.fh.fetch(reference=bait_chrom,start=left,end=right):
            items = item.split()
            ochrom, opos = items[2], int(items[3])
            targets[ochrom].append(opos)

        # count number of links
        counts = [0] *6
        self_cnt = 0
        mid = (left+right)/2
        for pos in targets[bait_chrom]:
            pos -= mid
            if   pos < -1000000:
                counts[0] += 1
            elif pos < -100000:
                counts[1] += 1
            elif pos < left-mid:
                counts[2] += 1
            elif pos < right-mid:
                self_cnt += 1
            elif pos < 100000:
                counts[3] += 1
            elif pos < 1000000:
                counts[4] += 1
            else:
                counts[5] += 1
        intra_cnt = sum(counts)
        chroms = sorted(targets,key=Utils.naturalkeys) 
        bait_links = pandas.DataFrame({'cnt':[len(targets[chrom]) for chrom in chroms],'chrom':chroms})
        bait_links.loc[bait_links.chrom==bait_chrom,'cnt'] = intra_cnt
        inter_cnt = bait_links.loc[bait_links.chrom!=bait_chrom,'cnt'].sum()

        # reads statistics
        Utils.touchtime("Parse read mapping rates and qualities ...")
        mdf = pandas.DataFrame({'Unmapped':[226,501],'Low MAPQ':[1056,816],'High MAPQ':[8721,8683]},index=['R1','R2'])

        # plotting
        Utils.touchtime("Plotting ...")
        import seaborn as sns
        import matplotlib.pyplot as plt
        sns.set_context('poster')
        sns.set_style('ticks')
        fig = plt.figure(figsize=[10,10])
        ax1 = plt.subplot2grid((20,10),(0,0),rowspan=10,colspan=4)
        ax2 = plt.subplot2grid((20,10),(0,6),rowspan=10,colspan=4)
        ax3 = plt.subplot2grid((20,10),(10,0),rowspan=3,colspan=10)
        ax4 = plt.subplot2grid((20,10),(15,0),rowspan=4,colspan=10)
        
        # read barplot
        Plot.ReadsBarPlot(mdf,ax1)

        # pie chart
        Plot.BaitCountPie([self_cnt,intra_cnt,inter_cnt],['Self-','Intra-','Inter-'],ax2)

        # bait count plot
        Plot.BaitCountPlot(sdepth,counts,ax3)

        # bait barplot
        Plot.BaitBarPlot(bait_links,bait_chrom,ax4)
        
        # save figure
        Utils.touchtime("Saving figure to {0} ...".format(outfile))
        plt.savefig(outfile)
        plt.close(fig)
    def GetIntraChromLinks(self,outfile=None,nbins=11,nperm=1000,seed=1024): 
        '''
        Get all local links from the provided region. 
        Get counts from the bait region to flanking N bins. The left and right bins with the same \
        distance to the bait region are merged.
        Parameters:
            outfile: string or None
                Save the count of intra-chromosome links.
            nbins: int
                Number of bins. The binsize is same as the bait peak size.
            nperm: int
                Number of permutations.
            seed: int
                random seed.
        '''
        chrom, size = self.bait_chrom, self.sizes[self.chroms.index(self.bait_chrom)]
        peaksize = self.peaksize
        left, right = self.left-peaksize/2, self.right-peaksize/2
        
        # permutation
        counts = []
        binsize = nbins*peaksize
        rs = numpy.random.RandomState(seed=seed)
        pcnt = 0
        starts = []
        while pcnt < nperm:
            start = rs.randint(0,size-binsize)
            end   = start + binsize
            mid = (start+end)/2
            if end < left or start > right:
                count = numpy.zeros(nbins,dtype=numpy.uint16)
                for item in self.fh.fetch(reference=chrom,start=start,end=end):
                    items = item.split()
                    pos, ochrom, opos = int(items[1]), items[2], int(items[3])
                    # check intra-chrom interactions
                    if ochrom==chrom and not start<opos<end: # same chrom
                        idx = min(abs(opos-start), abs(opos-end))/peaksize
                        if idx < nbins:
                            count[idx] += 1
                starts.append(start)
                counts.append(count)
                pcnt += 1
                if pcnt %1000 == 0:
                    Utils.touchtime("{0} permutations processed ...".format(pcnt))
        if nperm%1000:
            Utils.touchtime("{0} permutations processed ...     ".format(nperm))
        cdf = pandas.DataFrame(counts,columns=["bin_{0}".format(i+1) for i in range(nbins)])
        ns, ps = zip(*cdf.apply(Algorithms.NBFit,axis=0))
        if outfile:
            cdf['starts'] = starts
            cdf.to_csv(outfile,sep='\t',index=None)
        return ns, ps
    def GetInterChromLinks(self,outfile=None,binsize=1000000,nperm=1000,seed=1024):
        '''
        Get links between two chromosomal regions.
        Parameters:
            outfile: string
                file to save the permutation results
            binszie: int
                binsize to count links
            nperm: int
                number of permutations
            seed: int 
                seed
        Returns:
            inter_counts: numpy.array
                counts of inter-chrom links.
        '''
        inter_counts = numpy.zeros(nperm,dtype=numpy.uint16)
        chroms, sizes = zip(*[(chrom,size) for chrom,size in zip(self.chroms,self.sizes) if chrom!=self.bait_chrom and size>binsize])
        cumsizes = (numpy.array(sizes)-binsize).cumsum()
        bait_size = self.sizes[self.chroms.index(self.bait_chrom)]

        from bisect import bisect_left
        rs = numpy.random.RandomState(seed=seed)
        pcnt = 0
        while pcnt < nperm:
            # fetch bait_chrom random locus
            start = rs.randint(0,bait_size-binsize)
            end   = start + binsize
            if end < self.left or start > self.right:
                # fetch other chrom random locus
                pos = rs.randint(0,cumsizes[-1])
                idx = bisect_left(cumsizes,pos)
                ochrom, ostart, oend = chroms[idx], cumsizes[idx]-pos, cumsizes[idx]-pos+binsize
                for item in self.fh.fetch(reference=self.bait_chrom,start=start,end=end):
                    items = item.split()
                    pchrom, ppos = items[2], int(items[3])
                    # check inter-chrom interactions                
                    if ochrom==pchrom and ostart <= ppos < oend:
                        inter_counts[pcnt] += 1
                pcnt += 1
                if pcnt %1000 == 0:
                    Utils.touchtime("{0} permutations processed ...".format(pcnt))
        if nperm%1000:
            Utils.touchtime("{0} permutations processed ...     ".format(nperm))
        if outfile:
            numpy.savetxt(outfile,inter_counts,fmt="%d")
        # NB fit
        n, p = Algorithms.NBFit(inter_counts)
        return n, p

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
                    reads.add("{0}:{1}".format(items[6],items[7])) # read pair with the same starts only count once
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
            pos = int(pos)
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
        tbffile = "{0}.pairs".format(prefix)
        if overwrite or not (os.path.isfile(tbffile) or os.path.isfile(tbffile+".gz")):
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
        else:
            Utils.touchtime("Output file exists: {0}.pairs. Skipped ...".format(prefix))
        # creat tabix index
        if overwrite or not os.path.isfile(tbffile+".gz"):
            Utils.touchtime("Build index for {0}.pairs".format(prefix))
            tbffile = pysam.tabix_index("{0}.pairs".format(prefix),seq_col=0,start_col=1,end_col=1)
        else:
            tbffile += ".gz"
            Utils.touchtime("Index exists. Skipped ...")
        Utils.touchtime("FixMatePairs finished ...")
        return tbffile
    FixMatePairs=staticmethod(FixMatePairs)
    def DeterminePeakSize(depth,smooth_window=100):
        '''
        Determine the peak size of the bait region. The binsize is defined as the continous region around the bait (middle) position whose average window depth is larger than the average depth of the whole depth array.
        Paramegters:
            depth: numpy.array
                read depth
            smooth_window: int
                smooth window to determine the boundary.
        Returns:
            left, right: int
                boundary of the bait region.
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
        return left, right
    DeterminePeakSize=staticmethod(DeterminePeakSize)
    def RandomGenomicLoci(chroms,sizes,binsize=1000000, nloci=1000,seed=1024):
        '''
        Randomly sample positions from the whole genome.
        Parameters:
            chroms: list of strings
                chromosome names
            sizes: list of integers
                chromosome lengths
            binsize: integer
                preset binsize. Default is 1000,000.
            nloci: integer
                number of random genomic loci
            seed: integer
                Random seed used to initialize the pseudo-random number generator. 
        Yields:
            chrom: string
                chromosome name
            pos: integer
                genomic locus within [binsize/2, length-binsize/2]
        '''
        from bisect import bisect_left
        rs = numpy.random.RandomState(seed=seed)
        # remove short chromosomes
        schroms, ssizes = zip(*[(chrom,size) for chrom, size in zip(chroms,sizes) if size >binsize])
        cumsizes = (numpy.array(ssizes)-binsize).cumsum()
        for pos in rs.randint(0,cumsizes[-1],nloci):
            idx = bisect_left(cumsizes,pos)
            yield schroms[idx], cumsizes[idx]-pos+binsize/2
    RandomGenomicLoci=staticmethod(RandomGenomicLoci)
    def ReadStats(logfile):
        '''
        Parse read mapping rates and mapping qualities.
        '''
        with open(logfile) as fh:
            start = fh.tell()
            fh.seek(0,2) 
            fh.seek(max(start,fh.tell()-400)) 
            lines = fh.readlines()[-6:]
            total_reads = int(lines[0].split()[0])
            mapped_reads = total_reads - int(lines[2].split()[0])
        return total_reads, mapped_reads
    ReadStats=staticmethod(ReadStats)
    def NBFit(obs):
        '''
        '''
        df = pandas.Series(obs).value_counts()
        y, X = list(df.values), list(df.index)
        # fit by NB model
        mod = NBin(y,X)
        res = mod.fit()
        p, n = res.params
        return n, p
    NBFit=staticmethod(NBFit)

class Plot(object):
    def ReadsBarPlot(mdf,ax):
        '''
        '''
        mdf.plot.bar(stacked=True,ax=ax,legend=False,width=0.8,color=['whitesmoke','lightgrey','grey'])
        ax.legend(loc=2, bbox_to_anchor=(0.85, 0.9))
        ax.spines['right'].set_color('white')
        ax.spines['top'].set_color('white')
        ax.set_xticklabels(ax.get_xticklabels(),rotation=0)
        ax.set_ylabel("# of reads (M)")
        pdf = mdf.transpose()/mdf.sum(axis=1)*100
        y1,y2 = 0, 0
        for s1,s2,p1,p2 in zip(mdf.loc['R1',:],mdf.loc['R2',:],pdf['R1'],pdf['R2']):
            ax.text(0,y1+s1/2.,"{0:.1f}%".format(p1),ha='center')#,va='center')
            ax.text(1,y2+s2/2.,"{0:.1f}%".format(p2),ha='center')#,va='center')
            y1 += s1
            y2 += s2
        return ax
    ReadsBarPlot=staticmethod(ReadsBarPlot)
    def BaitCountPie(cnts,labels,ax):
        '''
        Pie chart of self-, intra- and inter-links.
        '''
        import matplotlib.pyplot as plt
        def make_autopct(values):
            def my_autopct(pct):
                total = sum(values)
                val = int(round(pct*total/100.0))
                return '{p:.2f}% ({v:d})'.format(p=pct,v=val)
            return my_autopct
        patches, texts, _ = ax.pie(cnts, labels=labels,colors= ['whitesmoke', 'lightblue', 'lightskyblue'],autopct=make_autopct(cnts), shadow=False, startangle=60, labeldistance=0.8,pctdistance=0.45,textprops={'ha':'left'})
        centre_circle = plt.Circle((0,0),0.75,color='black', fc='white',linewidth=1.25)
        ax.add_artist(centre_circle)
        ax.axis('equal')
        #ax.set_ylim(-1,1)
        #ax.set_xlim(-1,1)
        #ax.legend(patches,labels,ncol=3,loc='upper center')
        return ax
    BaitCountPie=staticmethod(BaitCountPie)
    def BaitCountPlot(sdepth,counts,ax):
        peaksize = sdepth.shape[0]
        mid = peaksize/2
        ax.plot(numpy.arange(-mid,-mid+peaksize)/float(peaksize),sdepth,linewidth=1)
        ax.set_xlim(-6.5,6.5)
        ax.spines['top'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.spines['right'].set_color('white')

        ax.set_yticks([])
        tick_loci = [-4.5,-2.5,-0.5,0,0.5,2.5,4.5]
        ax.set_xticks(tick_loci)
        ax.set_xticklabels(['-1M','-100K','','HS3','','100K','1M'])

        # counts
        counts_loci = [-5.5,-3.5,-1.5,1.5,3.5,5.5]
        ymin, ymax = ax.get_ylim()
        for c,i in zip(counts,counts_loci):
            ax.text(i,0.2*ymax,c,va='bottom')
        return ax
    BaitCountPlot=staticmethod(BaitCountPlot)
    def BaitBarPlot(df,bait_chrom,ax):
        '''
        Draw bait links barplot.
        '''
        # count of links from the bait region
        bait_cnt = df.cnt[df.chrom==bait_chrom].values[0]
        max_cnt = df.cnt[df.chrom!=bait_chrom].max()
        df.loc[df.chrom==bait_chrom,'cnt'] = int(max_cnt*1.5)

        # barplot
        ax.bar(df.index[::2],df.cnt[::2],width=0.9,color='darkgrey')
        ax.bar(df.index[1::2],df.cnt[1::2],width=0.9,color='lightgrey')
        ax.spines['top'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.spines['right'].set_color('white')
        ax.set_xlim(-0.5,df.index[-1]+0.5)
        ax.axis('off')

        # bait region
        ax.bar(df.index[df.chrom==bait_chrom],max_cnt*1.5,width=0.9,color='lightblue')
        ax.bar(df.index[df.chrom==bait_chrom],max_cnt*1.3,width=0.9,color='lightgrey')
        ax.bar(df.index[df.chrom==bait_chrom],max_cnt*1.2,width=0.9,color='lightblue')

        # numbers and chroms
        ymin, ymax = ax.get_ylim()
        offset = 0.03*(ymax-ymin)
        for i,chrom,cnt in zip(df.index,df.chrom,df.cnt):
            if chrom == bait_chrom:
                ax.text(i,cnt+offset,bait_cnt,ha='center',va='bottom')
            else:
                ax.text(i,cnt+offset,cnt,ha='center',va='bottom')
            ax.text(i,-offset,chrom[3:],ha='center',va='top')
        return ax
    BaitBarPlot=staticmethod(BaitBarPlot)

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
        return os.path.abspath(wdir)+"/"
    touchdir=staticmethod(touchdir)
    def naturalkeys(text):
        '''
        nlist = sorted(old_list,key=natural_keys) #sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        import re
        def atoi(text):
            return int(text) if text.isdigit() else text
        return [ atoi(c) for c in re.split('(\d+)', text) ]
    naturalkeys=staticmethod(naturalkeys)

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

class NBin(GenericLikelihoodModel):
    '''
    Estimate negative binomial parameters. Codes are modified from:
    http://www.statsmodels.org/dev/examples/notebooks/generated/generic_mle.html
    Usage:
        # generate random observations
        n, p = 5, 0.3
        obs = nbinom.rvs(n,p,size=1000)
        df = pandas.Series(obs).value_counts().sort_index()
        y, X = list(df.values), list(df.index)
        # fit by NB model
        mod = NBin(y,X)
        res = mod.fit()
        p, n = res.params
    '''
    def __init__(self, endog, exog, **kwds):
        super(NBin, self).__init__(endog, exog, **kwds)
    def _ll_nb2(y, X, beta, alph):
        mu = numpy.exp(numpy.dot(X, beta))
        size = 1/alph
        prob = size/(size+mu)
        ll = nbinom.logpmf(y, size, prob)
        return ll
    _ll_nb2=staticmethod(_ll_nb2)        
    def nloglikeobs(self, params):
        alph = params[-1]
        beta = params[:-1]
        ll = NBin._ll_nb2(self.endog, self.exog, beta, alph)
        return -ll     
    def fit(self, start_params=None, maxiter=10000, maxfun=5000, disp=0, **kwds):
        # we have one additional parameter and we need to add it for summary
        self.exog_names.append('alpha')
        if start_params == None:
            # Reasonable starting values
            start_params = numpy.append(numpy.zeros(self.exog.shape[1]), .5)
            # intercept
            start_params[-2] = numpy.log(self.endog.mean())
        return super(NBin, self).fit(start_params=start_params, 
                                     maxiter=maxiter, maxfun=maxfun, 
                                     disp=disp, **kwds) 

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example:"+sys.argv[0]+" genome fqfile prefix ")
    Tools.bowtie2_SE(sys.argv[1], sys.argv[2],sys.argv[3])

