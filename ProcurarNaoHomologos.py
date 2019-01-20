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
        """Módulo que recebe como parâmetros um ficheiro com formato xml e um ficheiro de output formato txt"""
        handle=open(ficheiro,"r")
        blast_records=NCBIXML.parse(handle)
        evalue=0.05
        for blast_record in blast_records:
            algn=[]  #lista com alinhamentos
            for alignment in blast_record.alignments:
                res=0
                for hsp in alignment.hsps:
                    if hsp.expect>evalue:
                        res+=1
                if res==len(alignment.hsps):  #se todos os high scoring pairs forem inferiores a 0.05, então o alinhamento é adicionado
                    algn.append(blast_record.query)
            if len(algn)==len(blast_record.alignments):  #se todos os alinhamentos tiverem high scoring pairs inferiores a 0.05, então adiciona-se a proteína, pois não é homóloga
                self.__nao_homologos.append(blast_record.query)
        fich=open(out,"w")
        for gene in self.__nao_homologos:
            fich.write(gene+"\n")
        fich.close()
        handle.close()
    
    def print_nao_homologos(self):
        print(" ".join(self.__nao_homologos))
        