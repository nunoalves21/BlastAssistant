# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 12:00:33 2019

@author: mario
"""
from Bio.Blast import NCBIXML

class ProcurarNaoHomologos:
    def __init__(self):
        self.__nao_homologos=[]
    
    def procura_naohomologos(self,ficheiro,out):
        handle=open(ficheiro,"r")
        blast_records=NCBIXML.parse(handle)
        evalue=0.05
        for blast_record in blast_records:
            algn=[]
            for alignment in blast_record.alignments:
                res=0
                for hsp in alignment.hsps:
                    if hsp.expect>evalue:
                        res+=1
                if res==len(alignment.hsps):
                    algn.append(blast_record.query)
            if len(algn)==len(blast_record.alignments):
                self.__nao_homologos.append(blast_record.query)
        fich=open(out,"w")
        for gene in self.__nao_homologos:
            fich.write(gene+"\n")
        fich.close()
        handle.close()
    
    def print_nao_homologos(self):
        print(" ".join(self.__nao_homologos))