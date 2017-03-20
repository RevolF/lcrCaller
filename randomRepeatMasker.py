# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:41:30 2017

@author: thor
"""

import os
from optparse import OptionParser
import subprocess as sp
import random


parser=OptionParser()
parser.add_option(
    '-T',
    '--work-dir',
    dest='workDir',
    help='work directory'
    )
    
parser.add_option(
    '-R',
    '--ref-gene',
    dest='refGene',
    help='ref gene file, default set to /home/ljzhang/data/refGene/refGene.txt',
    default='/home/ljzhang/data/refGene/refGene.txt'
    )
    
parser.add_option(
    '-N',
    '--random-nbr',
    dest='randNbr',
    help='number of gene for random choices'
    )
    
parser.add_option(
    '-O',
    '--output-file',
    dest='outFile',
    help='output file name'
    )
    
parser.add_option(
    '-P',
    '--cpu-nbr',
    dest='cpuNbr',
    help='cpu nbrs, default set to 1',
    default='1'
    )
    
(options,args)=parser.parse_args()
if not options.workDir or not options.refGene or not options.randNbr or not options.outFile:
    parser.print_help()
    exit(1)
    
assert os.path.exists(options.workDir) and os.path.exists(options.refGene)

def main():
    global options
    os.chdir(options.workDir)
    idxToPosLst,candRefGeneCount=getRandomGeneInfo(options)
    assert int(options.cpuNbr) < candRefGeneCount
    for idxItem in idxToPosLst:
        runRepeatMasker(options.workDir,idxItem[0],idxItem[1],idxItem[2],idxItem[3],options.cpuNbr)
    os.chdir(options.workDir)
    return
    
def getRandomGeneInfo(options):
    lineCount=0
    if not os.path.exists(options.workDir+'/refGene_20moreExons.txt'):
        awkCmd="awk '$9>20{print $0}' %s > %s" % (options.refGene,options.workDir+'/refGene_20moreExons.txt')
        sp.call(awkCmd,shell=True)
    with open(options.workDir+'/refGene_20moreExons.txt','r') as infh:
        for line in infh.xreadlines():
            lineCount+=1
            pass
    
    print lineCount
    randIdxDct={}
    for i in range(int(options.randNbr)):
        randIdx=random.randint(0,lineCount)
        while randIdxDct.has_key(randIdx):
            randIdx=random.randint(0,lineCount)
        randIdxDct[randIdx]=1
    
    randIdxLst=sorted(randIdxDct.keys())
    print randIdxLst
	
    idxToPosLst=[]
    ind=0
    with open(options.workDir+'/refGene_20moreExons.txt','r') as infh, open(options.workDir+'/randRefGenes.txt','w') as outfh:
        for line in infh.xreadlines():
            if ind in randIdxLst:
                print ind
                linear=line.strip().split('\t')
                idxToPosLst.append([linear[1],linear[2],linear[4],linear[5]])
                outfh.write(line)
            ind+=1
    return idxToPosLst,lineCount
    
def runRepeatMasker(workDir,trxName,chr,startPos,endPos,cpuNbr):
    if not os.path.exists(workDir+'/'+trxName):
        os.mkdir(workDir+'/'+trxName)
    os.chdir(workDir+'/'+trxName)
    samfaidxCmd='samtools faidx /media/disk2/ljzhang/data/hg19/ucsc.hg19.fasta %s > %s' % (chr+':'+startPos+'-'+endPos,trxName+'.fasta')
    sp.call(samfaidxCmd,shell=True)
#    ~/bin/RepeatMasker/RepeatMasker -species human dmdRegion_10k.fasta -a -inv -lcambig -xsmall -poly -source -html -ace -gff -u -xm -pa 12
    repMskrCmd='/home/ljzhang/bin/RepeatMasker/RepeatMasker -species human %s -a -inv -lcambig -xsmall -poly -source -html -ace -gff -u -xm -pa %s' % (trxName+'.fasta',cpuNbr)
    sp.call(repMskrCmd,shell=True)
    return

def lcrFinder():
    
    pass

if __name__=='__main__':
    main()



