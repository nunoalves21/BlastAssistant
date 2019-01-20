# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 15:28:47 2019

"""
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo

class AlignmentReader:
    def __init__(self,file):
        self.__alignment=AlignIO.read(file,"fasta")
        
    def getAlignment(self):
        return self.__alignment
    
    def getConsenso(self):
        align=SummaryInfo(self.__alignment)
        return align.gap_consensus()
    
    def getConservedDomain(self):
        cons=[]
        align=SummaryInfo(self.__alignment)
        consenso=str(align.gap_consensus())
        temp=''
        for i in range(len(consenso)):
            if consenso[i] not in "X-":
                temp+=consenso[i]
            else:
                if temp!='':
                    cons.append(temp)
                temp=''
        max_cons=''
        for i in cons:
            if len(i)>len(max_cons):max_cons=i
        return max_cons
        
    def findConservedDomain(self):
        return self.getConsenso().find(self.getConservedDomain())
    
    def getPerc(self):
        align=SummaryInfo(self.__alignment)
        return float(1- (align.gap_consensus().count("X")+align.gap_consensus().count("-"))/len(str(align.gap_consensus())))
        
        
        
        
        
        
        