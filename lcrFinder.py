# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 12:22:10 2017

@author: thor
"""

# this script is used for searching for potential LCR structure based on detecting result of RepeatMasker

from __future__ import division
import os,glob
import subprocess as sp
from optparse import OptionParser

parser=OptionParser()

parser.add_option(
	'-D',
	'--raw-dir',
	dest='destDir',
	help='destDir: /media/disk2/ljzhang/project/dmd/potentialLcrFinder',
	default='/media/disk2/ljzhang/project/dmd/potentialLcrFinder'
	)
	
parser.add_option(
	'-F',
	'--csv-file',
	dest='csvFile',
	help='dmdRegion_10k.fasta.out.exonAnno.sorted.4456.exon.spanning.sorted.indexed.csv',
	default='dmdRegion_10k.fasta.out.exonAnno.sorted.4456.exon.spanning.sorted.indexed.csv'
	)
	
parser.add_option(
	'-W',
	'--win-size',
	dest='winSize',
	help='windows size , default set to 15',
	default='15'
	)
	
(options,args)=parser.parse_args()

destDir=options.destDir
csvFile=options.csvFile
winSize=int(options.winSize)

#==============================================================================
# destDir='/media/disk2/ljzhang/project/dmd/potentialLcrFinder'
# csvFile='dmdRegion_10k.fasta.out.exonAnno.sorted.4456.exon.spanning.sorted.indexed.csv'
# winSize=15
#==============================================================================

def main():
	global destDir,csvFile,winSize
	os.chdir(destDir)
	idx2Indt,idx2AllInfo=getDctReady(destDir,csvFile)
	'''
	idx2Indt:
		'1' => 'ALUY'
	idx2AllInfo:
		'1' => ['ALUY','31528625','31528926','C','SINE/Alu','INTRON56_55']
		
	sortedIdxKeys:
		['1','2',...]
		[identifierIdx keys for splitting]
	'''
	sortedIdxKeys=sorted(idx2Indt.keys(),key=lambda x:int(x))
	
	if not os.path.exists(destDir+'/'+str(winSize)):
		os.mkdir(destDir+'/'+str(winSize))
	
	os.chdir(destDir+'/'+str(winSize))
	
	for slide in range(0,winSize):
		mainExc(winSize,sortedIdxKeys,slide,idx2Indt,destDir)
	
	resFiles=glob.glob('*results')
	outfile='potentialLcrMerged.txt'
	for resFile in resFiles:
		sp.call('cat '+resFile+' >> '+outfile,shell=True)
	
	return
	
def mainExc(winSizeMe,sortedIdxKeysMe,slideMe,idx2IndtMe,destDirMe):
	
	os.chdir(destDirMe+'/'+str(winSizeMe))
	outfh=open(str(slideMe)+'.results','w')
	
	winKeyDct=split2WinSize(winSizeMe,sortedIdxKeysMe,slideMe)
	winKeyDctKeys=sorted(winKeyDct.keys(),key=lambda x:int(x))
	'''
	winKeyDct:
		'0' => ['1','2',...,'15']
		windowIdx => [identifierIdx]
	winKeyDctKeys:
		['0','1','2','3']
		windowIdx
	'''
	cmpedKeySet=set()
	for i in winKeyDctKeys:
		cmpedKeySet.add(i)
		for j in winKeyDctKeys:
			if j in cmpedKeySet:
				continue
			score, primeMatIdxLst, subMatIdxLst = interCmp(winKeyDct,i,j,idx2IndtMe)
			primeMatIdxLst = sorted(primeMatIdxLst,key=lambda x: int(x))
			subMatIdxLst = sorted(subMatIdxLst,key=lambda x: int(x))
			'''
			primeMatIdxLst, subMatIdxLst:
				identifiers
			'''
			outfh.write('win_'+str(winSizeMe)+'\tslide_'+str(slideMe)+'\t'+i+'_'+j+'\t'+str(score)+'\t'+str(round(score/winSizeMe,4))+'\tPrime\t')
			outfh.write('_'.join(primeMatIdxLst)+'\tEnd\t')
			outfh.write('_'.join(subMatIdxLst)+'\tPrimeIdentifier\t')
			outfh.write('_'.join([idx2IndtMe[ai] for ai in primeMatIdxLst])+'\tEndIdentifier\t')
			outfh.write('_'.join([idx2IndtMe[aj] for aj in subMatIdxLst])+'\n')
			
	outfh.close()
	return

def getDctReady(destDirGd,csvGd):
	'''
	idx2Indt:
		'1' => ALUY
	idx2AllInfo:
		'1' => ['ALUY','31528625','31528926','C','SINE/Alu','INTRON56_55']
	'''
	os.chdir(destDirGd)
	infh=open(csvGd,'r')
	idx2Indt={}
	idx2AllInfo={}
	
	for line in infh.xreadlines():
		linear=line.strip().split(',')
		idx2Indt[linear[0]]=linear[1]
		idx2AllInfo[linear[0]]=linear[1:]
		
	return idx2Indt,idx2AllInfo
	
def split2WinSize(winSizeSw,sortedIdxKeysSw,slideSw):
	sortedIdxKeysSw=sortedIdxKeysSw[slideSw:]
	keyLengh=len(sortedIdxKeysSw)
	compleSize=keyLengh//winSizeSw
	
	winKeyDct={}
	end=0
	for i in range(compleSize):
		stt=i*winSizeSw
		end=stt+winSizeSw
		winKeyDct[str(i)]=sortedIdxKeysSw[stt:end]
	if end < len(sortedIdxKeysSw):
		winKeyDct[str(compleSize)]=sortedIdxKeysSw[end:len(sortedIdxKeysSw)]
	
	return winKeyDct
	
def interCmp(winKeyDctIc,primeIc,subIc,idx2IndtIc):
	'''
	winKeyDct:
		'0' => ['1','2',...,'15']
		windowIdx => [identifierIdx]
	primeLst, subLst;
		identifier list
	'''
	primeLst=winKeyDctIc[primeIc]
	subLst=winKeyDctIc[subIc]
	score=0
	primeMatchedIdxList=[]
	subMatchedIdxList=[]

	for primeIdx in primeLst:
		primeItem=idx2IndtIc[primeIdx]
		cmpSet=set()
		for subIdx in subLst:
			subItem=idx2IndtIc[subIdx]
			if subItem in cmpSet:
				continue
			if primeItem == subItem:
				score += 1
				primeMatchedIdxList.append(primeIdx)
				subMatchedIdxList.append(subIdx)
				cmpSet.add(subItem)
	return score, primeMatchedIdxList, subMatchedIdxList
	
if __name__=='__main__':
	main()




