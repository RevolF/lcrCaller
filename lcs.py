# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 17:08:32 2017

@author: thor
"""

# implementation of LCS algos regarding to longest common subsequence

class Lcs:
    def __init__(self,longStr,shortStr):
        '''
        note that long string can be longer and equal to the size of short string
        in returning results:
            self.lcs : longest common sequence list
            self.refLong : returns refined long string in list form, '_' stands for deletion
            self.refShort : returns refined short string in list form, '_' stands for deletion
        '''
        self.longStr=longStr
        self.shortStr=shortStr
        self.colNbr=len(longStr)
        self.rowNbr=len(shortStr)
        self._matrix=''
        self.lcs=[]
        self.refLong=[]
        self.refShort=[]
        
        
    def _initializeMatrix(self):
        '''
        initialize matrix: each element of longStr as index of columns
                and shortStr for index of rows
        '''
        cellMatrix=[]
        firstRow=[Cell('NASTR','NASTR',-1,-1)]
        for i in range(self.colNbr):
            firstRow.append(Cell('NASTR',self.longStr[i],0,i))
        cellMatrix.append(firstRow)
        for i in range(self.rowNbr):
            tmpRow=[Cell(self.shortStr[i],'NASTR',i,-1)]
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
                    self.lcs.append(self._matrix[i][j].lowerEle)
                else:
                    self._matrix[i][j].value = min(self._matrix[i-1][j],self._matrix[i-1][j-1],self._matrix[i][j-1]) + 1
        return
    
    def _backTrace(self):
        i=self.rowNbr
        j=self.colNbr
        curCell=self._matrix[i][j]
        
        self.refLong=[self._matrix[i][j].upperEle]
        self.refShort=[self._matrix[i][j].lowerEle]
        
        while i>=0 and j>=0:
            cellPointer = self._returnPointer(self._matrix[i-1][j-1],self._matrix[i-1][j],self._matrix[i][j-1])
            if cellPointer == 1:
                curCell.prevCell=self._matrix[i-1][j-1]
                self.refLong.insert(0,self._matrix[i-1][j-1].upperEle)
                self.refShort.insert(0,self._matrix[i-1][j-1].lowerEle)
                i -= 1
                j -= 1
            elif cellPointer == 2:
                curCell.prevCell=self._matrix[i-1][j]
                self.refLong.insert(0,'_')
                self.refShort.insert(0,self._matrix[i-1][j].lowerEle)
                i -= 1
            elif cellPointer == 3:
                curCell.prevCell=self._matrix[i][j-1]
                self.refLong.insert(0,self._matrix[i][j-1].upperEle)
                self.refShort.insert(0,'_')
                j -= 1
            curCell=curCell.prevCell
            
        return
    
    def _returnPointer(t1,t2,t3):
        '''
        note that t1,t2,t3 are all Cell objects
            t1 indicates topleft
            t2 indicates topright
            t3 indicates bottomleft
        for cellPointer:
            1 indicates topleft
            2 indicates topright
            3 indicates bottomleft
        '''
        if t1.value <= t2.value and t1.value <= t3.value:
            return 1
        elif t2.value < t1.value and t2.value < t3.value:
            return 2
        elif t3.value < t1.value and t3.value <= t2.value:
            return 3
        
    def analysis(self):
        self._initializeMatrix()
        self._processMatrix()
        self._backTrace()
        del self._matrix
        return
        


class Cell:
    def __init__(self,lowerEle,upperEle,rowNbr,colNbr,value=0,cellPointer='',prevCell=''):
        '''
        for cellPointer:
            1 indicates topleft
            2 indicates topright
            3 indicates bottomleft
        '''
        self.upperEle=upperEle
        self.lowerEle=lowerEle
        self.rowNbr=rowNbr
        self.colNbr=colNbr
        self.value=value
        self.cellPointer=cellPointer
        self.prevCell=prevCell
    
    
