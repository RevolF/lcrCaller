# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:28:00 2017

@author: thor
"""

# implementation of lcr_finder using lcs script

from __future__ import division
import re,sys
from optparse import OptionParser

parser=OptionParser()

parser.add_option(
    '-D',
    '--work-dir',
    dest='workDir',
    help='working dir, all results file will be put here'
    )
    
parser.add_option(
    '-R',
    '--rpout-file',
    dest='rpFile',
    help='repeat masker output file'
    )
    
'''
    put lcs class inside the script
'''
#==============================================================================
# parser.add_option(
#     '-L',
#     '--lcs-script',
#     dest='lcsScpt',
#     help='full path to lcs.py',
#     default='/home/ljzhang/project/lcs.py'
#     )
#==============================================================================
    
parser.add_option(
    '-E',
    '--exon-dict',
    dest='exonDct',
    help='exon dict for exon annotation, default set to /media/disk2/ljzhang/data/dmdExons/dmd.exons.txt',
    default='/media/disk2/ljzhang/data/dmdExons/dmd.exons.txt'
    )

parser.add_option(
    '-W',
    '--win-size',
    dest='winSize',
    help='windows size for sliding, by default this is set to 20',
    default='20'
    )
    
(options,args)=parser.parse_args()

if not options.workDir or not options.rpFile:
    parser.print_help()
    sys.exit(1)

#==============================================================================
# if not os.path.exists(options.lcsScpt):
#     print 'lcs.py(c) not found in the given parameters'
#     parser.print_help()
#     sys.exit(1)
#==============================================================================

def main():
    global options
    exonDct=getExonDct(options.exonDct)
    infoMat=splitRpFile(options.rpFile,exonDct)
    
    rpName=options.rpFile.split('/')[-1].rstrip('out')
    
    matOutName=options.workDir + '/' + rpName + 'mat.out'
    writeMatrix(infoMat,matOutName)
    
    lcrResFile=options.workDir+'/' + rpName + 'winSize.' + options.winSize + '.lcsResMat.out'
    slideWin(infoMat,options.winSize,lcrResFile,exonDct)
    return
    
def splitRpFile(rpFile, exonDct):
    '''
    splitting repeat masker .out file and gives:
        start position on chrX
        end position on chrX
        matching repeat
        repeat class/family
        ID
    using cutoff:
        perc div. : 15%
        perc del. : 8%
        perc ins. : 8%
    '''
    infh=open(rpFile,'r')
    infoMat=[]
    
    ind = 0
    for line in infh.xreadlines():
        line=line.strip()
        if ind < 3:
            ind += 1
            continue
        linear=re.split('\s+',line)
        '''
        magic variables here, change if plan modified
        '''
        if float(linear[1]) < 15 and float(linear[2]) < 8 and float(linear[3]) < 8:
            markStt=linear[4].split('-')[0].split(':')[1]
            sttPos=int(markStt)+int(linear[5])
            endPos=int(markStt)+int(linear[6])
            orient=linear[8]
            eleId=linear[9]
            eleFam=linear[10]
            lineId=linear[14]
            regionAnno=exonAnno(exonDct, sttPos)
            infoMat.append([sttPos,endPos,orient,eleId,eleFam,lineId,regionAnno])
    return infoMat

def getExonDct(dmdExon):
    exonDct={}
    infh=open(dmdExon,'r')
    intronLst=[]
    intronName='INTRON79'
    ind=1
    for line in infh.xreadlines():
        linear=line.strip().split('\t')
        rawName=linear[0].lstrip('EXON')
        if ind == 1:
            intronLst.append(linear[2])
            exonDct[linear[0]]=[linear[1],linear[2]]
            ind += 1
            continue
        intronLst.append(linear[1])
        intronName = intronName + '_' + rawName
        exonDct[intronName] = intronLst
        intronName = 'INTRON' + rawName
        intronLst = [linear[2]]
        exonDct[linear[0]] = [linear[1],linear[2]]
    infh.close()
    return exonDct

def writeMatrix(infoMat,fileName):
    outfh=open(fileName,'w')
    for subList in infoMat:
        outfh.write('\t'.join([str(i) for i in subList])+'\n')
    outfh.close()
    return

def slideWin(infoMat,winSize,lcsResFile,exonDct):
    winSize=int(winSize)
    
    print 'inside slideWin'
    
    outfh=open(lcsResFile,'w')
    header='winStartPos\twinFirstId\twinNextId\tmatchScore\trepEleLcs\twinFirstLcsOrient\twinNextLcsOrient\twinFirstString\twinNextString\twinFirstRegion\twinNextRegion\n'
    outfh.write(header)
    
    rowNbr=len(infoMat)
    
    for i in range(winSize):
        repEleMat=returnRepEleMat(infoMat,rowNbr,i,winSize)
        print 'winSlide: '+str(i)+' out of '+str(winSize)
        repEleMatRowNbr=len(repEleMat)
        
        for j in range(repEleMatRowNbr):
            if (j+1) == repEleMatRowNbr:
                break
            for k in range(j+1,repEleMatRowNbr):
                lcsObj=Lcs(repEleMat[j],repEleMat[k])
                score=round(200*len(lcsObj.lcs)/(lcsObj.colNbr+lcsObj.rowNbr),3)
                
                winFirstLcsOrient=''.join(lcsObj.upperLcsOri)
                winNextLcsOrient=''.join(lcsObj.lowerLcsOri)
                winFirstString=','.join(lcsObj.refLong)
                winNextString=','.join(lcsObj.refShort)
                
                winFirstRegion=winRegionAnno(exonDct,lcsObj.upperLcsSttPos)
                winNextRegion=winRegionAnno(exonDct,lcsObj.lowerLcsSttPos)
                
                writeList=[i,j,k,score,lcsObj.lcs,winFirstLcsOrient,winNextLcsOrient,winFirstString,winNextString,winFirstRegion,winNextRegion]
                
                outfh.write('\t'.join([str(subItem) for subItem in writeList])+'\n')
    outfh.close()
    return
    
def winRegionAnno(exonDct,sttPosLst):
    sttPosLst=[int(i) for i in sttPosLst]
    if sttPosLst==[]:
        return 'NASTR'
    sttPos=min(sttPosLst)
    endPos=max(sttPosLst)
    sttPosAnno=exonAnno(exonDct,sttPos)
    endPosAnno=exonAnno(exonDct,endPos)
    winRegion=sttPosAnno+'---'+endPosAnno
    return winRegion
  
def returnRepEleMat(infoMat,rowNbr,sttRow,winSize):
    '''
    repEleMat is a matrix containing rows of elements of RepeatEle
    '''
    repEleMat=[]
    curRow=sttRow
    while (curRow + winSize) <= rowNbr:
        tmpList=[]
        for i in range(curRow,curRow+winSize):
            '''
            class RepeatEle:
                def __init__(self,repName,repOri,repId,repRegion):
            [sttPos  endPos  orient  eleId  eleFam  lineId  regionAnno]
            '''
            tmpList.append(RepeatEle(infoMat[i][3],infoMat[i][2],infoMat[i][5],infoMat[i][0]))
        repEleMat.append(tmpList)
        curRow += winSize
    if curRow < rowNbr:
        tmpList=[]
        for i in range(curRow,rowNbr):
            tmpList.append(RepeatEle(infoMat[i][3],infoMat[i][2],infoMat[i][5],infoMat[i][0]))
        repEleMat.append(tmpList)
    return repEleMat
    

def exonAnno(exonDct,sttPos):
    for key in exonDct.keys():
        if int(exonDct[key][0]) < int(sttPos) and int(sttPos) <= int(exonDct[key][1]):
            return key
    return 'NASTR'




class Lcs:
    def __init__(self,longStr,shortStr):
        '''
        note that long string can be longer and equal to the size of short string
        in returning results:
            self.lcs : longest common sequence list
            self.refLong : returns refined long string in list form, '_' stands for deletion
            self.refShort : returns refined short string in list form, '_' stands for deletion
        adding a region start and end position to Lcs, by extracting max and min of coordinates
        '''
        if len(longStr) >= len(shortStr):
            self.longStr=longStr
            self.shortStr=shortStr
        else:
            self.longStr=shortStr
            self.shortStr=longStr
        
        self.colNbr=len(self.longStr)
        self.rowNbr=len(self.shortStr)
        self._matrix=''
        self.lcs=[]
        self.refLong=[]
        self.refShort=[]
        
        self.lcs=[]
        self.upperLcsOri=[]
        self.lowerLcsOri=[]
        self.upperLcsSttPos=[]
        self.lowerLcsSttPos=[]
        
        self._analysis()
        
    def _initializeMatrix(self):
        '''
        initialize matrix: each element of longStr as index of columns
                and shortStr for index of rows
        Cell:
            lowerEle,upperEle,rowNbr,colNbr,value=0,cellPointer='',prevCell=''
        RepeatEle:
            self,repName,repOri,repId,sttPos
        '''
        cellMatrix=[]
        emptyRepEle=RepeatEle('NASTR','NASTR','NASTR','0')
        firstRow=[Cell(emptyRepEle,emptyRepEle,-1,-1)]
        for i in range(self.colNbr):
            firstRow.append(Cell(emptyRepEle,self.longStr[i],0,i,i+1))
        cellMatrix.append(firstRow)
        for i in range(self.rowNbr):
            tmpRow=[Cell(self.shortStr[i],emptyRepEle,i,-1,i+1)]
            for j in range(self.colNbr):
                tmpRow.append(Cell(self.shortStr[i],self.longStr[j],i,j))
            cellMatrix.append(tmpRow)
        self._matrix=cellMatrix
        return
    
    def _processMatrix(self):
        for i in range(1,self.rowNbr+1):
            for j in range(1,self.colNbr+1):
                if self._matrix[i][j].lowerEle == self._matrix[i][j].upperEle:
                    self._matrix[i][j].value = self._matrix[i-1][j-1].value
                    '''
                    self.lcs.append(self._matrix[i][j].lowerEle)
                     ### note that this would cause false results inside a for loop ###
                        use a cellPointer for indicating lcs
                    '''
                    
                else:
                    self._matrix[i][j].value = min(self._matrix[i-1][j].value,self._matrix[i-1][j-1].value,self._matrix[i][j-1].value) + 1
        return
    
    def _backTrace(self):
        i=self.rowNbr
        j=self.colNbr
        curCell=self._matrix[i][j]
        
        '''
        when applying backtracing, should place string after condition returns true
        '''
        '''
        self.refLong=[self._matrix[i][j].upperEle]
        self.refShort=[self._matrix[i][j].lowerEle]
        '''
        
        self.refLong=[]
        self.refShort=[]
        '''
        avoid loop adding string when there's a blanc
        '''
        
        while i>=1 and j>=1:
            cellPointer = self._returnPointer(self._matrix[i-1][j-1],self._matrix[i-1][j],self._matrix[i][j-1])
            if cellPointer == 1:
                curCell.prevCell=self._matrix[i-1][j-1]
#                self.refLong.insert(0,self._matrix[i-1][j-1].upperEle)
#                self.refShort.insert(0,self._matrix[i-1][j-1].lowerEle)
                self.refLong.insert(0,self._matrix[i][j].upperEle)
                self.refShort.insert(0,self._matrix[i][j].lowerEle)
                if self._matrix[i][j].upperEle == self._matrix[i][j].lowerEle:
                    self.lcs.insert(0,self._matrix[i][j].upperEle)
                    
                    '''
                    Cell:
                        lowerEle,upperEle,rowNbr,colNbr,value=0,cellPointer='',prevCell=''
                    '''
                    
                    self.upperLcsOri.insert(0,self._matrix[i][j].upperOri)
                    self.lowerLcsOri.insert(0,self._matrix[i][j].lowerOri)
                    
                    self.upperLcsSttPos.insert(0,self._matrix[i][j].upperSttPos)
                    self.lowerLcsSttPos.insert(0,self._matrix[i][j].lowerSttPos)
                    
                i -= 1
                j -= 1
            elif cellPointer == 2:
                curCell.prevCell=self._matrix[i-1][j]
                self.refLong.insert(0,'_')
#                self.refShort.insert(0,self._matrix[i-1][j].lowerEle)
                self.refShort.insert(0,self._matrix[i][j].lowerEle)
                i -= 1
            elif cellPointer == 3:
                curCell.prevCell=self._matrix[i][j-1]
#                self.refLong.insert(0,self._matrix[i][j-1].upperEle)
                self.refLong.insert(0,self._matrix[i][j].upperEle)
                self.refShort.insert(0,'_')
                j -= 1
            curCell=curCell.prevCell
            
        if i == 0 and j >= 1:
            while j >= 1:
                self.refLong.insert(0,self._matrix[i][j].upperEle)
                self.refShort.insert(0,'_')
                j -= 1
        elif j == 0 and i >= 1:
            while i >= 1:
                self.refShort.insert(0,self._matrix[i][j].lowerEle)
                self.refLong.insert(0,'_')
                i -= 1
            
        return
    
    def _returnPointer(self,t1,t2,t3):
        '''
        note that t1,t2,t3 are all Cell objects
            t1 indicates topleft
            t2 indicates topright
            t3 indicates bottomleft
        for cellPointer:
            1 indicates topleft
            2 indicates topright
            3 indicates bottomleft
        should notice a special case:
            backtracing has reached the first row or col
        '''
        if t1.value <= t2.value and t1.value <= t3.value:
            return 1
        elif t2.value < t1.value and t2.value < t3.value:
            return 2
        elif t3.value < t1.value and t3.value <= t2.value:
            return 3
        
    def _analysis(self):
        self._initializeMatrix()
        self._processMatrix()
        self._backTrace()
#        del self._matrix
        return
        
    def getMatrix(self):
        matrix=[]
        for i in range(self.rowNbr+1):
            tmpRow=[a.value for a in self._matrix[i]]
            matrix.append(tmpRow)
        return matrix


class Cell:
    def __init__(self,lowerEle,upperEle,rowNbr,colNbr,value=0,cellPointer='',prevCell=''):
        '''
        for cellPointer:
            1 indicates topleft
            2 indicates topright
            3 indicates bottomleft
        add a flag in Cell to ensure that each element is added only once
        '''
        self.upperEle=upperEle.repName
        self.lowerEle=lowerEle.repName
        
        self.upperOri=upperEle.repOri
        self.lowerOri=lowerEle.repOri
        
        self.upperId=upperEle.repId
        self.lowerId=lowerEle.repId
        
        self.upperSttPos=upperEle.sttPos
        self.lowerSttPos=lowerEle.sttPos
        
        self.rowNbr=rowNbr
        self.colNbr=colNbr
        self.value=value
        self.cellPointer=cellPointer
        self.prevCell=prevCell

class RepeatEle:
    def __init__(self,repName,repOri,repId,sttPos):
        self.repName=repName
        self.repOri=repOri
        self.repId=repId
        self.sttPos=sttPos
        
        
if __name__=='__main__':
    main()
