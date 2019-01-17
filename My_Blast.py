# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 10:49:10 2019

@author: mario
"""
import time 
from Bio.Blast import NCBIWWW
from Bio import SeqIO

class My_Blast:
    def __init__(self,file_prot,file_out,db="nr"):
        self.__file_prot=file_prot
        self.__out=file_out
        self.__db=db
        
    def make_blast(self):
        records=SeqIO.parse(self.__file_prot,"fasta")
        save_file=open(self.__out,"w")
        for record in records:
            now=time.time()
            result_handle = NCBIWWW.qblast("blastp", self.__db, record.format("fasta"), entrez_query='Homo sapiens [organism]')
            save_file.write(result_handle.read()+"\n")
            end=time.time()
            print("A proteína %s já foi submetida ao blast e demorou %s segundos. "%(record.id,end-now))
        save_file.close()
        records.close()
        