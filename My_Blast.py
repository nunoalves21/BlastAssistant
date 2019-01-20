# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 10:49:10 2019

"""
import time 
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

class My_Blast:
    def __init__(self,file_prot,file_out,db="nr"):
        """O construtor recebe como parâmetros um ficheiro em formato fasta, um ficheiro de output em formato xml e a base de dados a utilizar"""
        self.__file_prot=file_prot
        self.__out=file_out
        self.__db=db
        
    def make_blast(self):
        """Faz um blast das proteínas que se encontram no ficheiro em formato fasta contra o genoma humano e imprime o tempo de duração"""
        records=SeqIO.parse(self.__file_prot,"fasta")
        save_file=open(self.__out,"w")
        for record in records:
            beginning=time.time()
            result_handle = NCBIWWW.qblast("blastp", self.__db, record.format("fasta"), entrez_query='Homo sapiens [organism]')
            save_file.write(result_handle.read()+"\n")
            end=time.time()
            print("A proteína %s já foi submetida ao blast e demorou %s segundos. "%(record.id,end-beginning))
        save_file.close()
        records.close()
    
    def make_blast_completar(self):
        """Módulo que procura a última proteína a ser submetida a um blast no ficheiro de formato xml e continua o blast a partir dessa."""
        result_handle=open(self.__out)
        blast_records=NCBIXML.parse(result_handle)
        blasts=[]
        for blast_record in blast_records: 
            blasts.append(blast_record.query)
        blast_records.close()
        last_prot=blasts[len(blasts)-1]
        records=SeqIO.parse(self.__file_prot,"fasta")
        save_file=open(self.__out,"a")
        save_file.write("\n")
        encontrada=False
        for record in records:
            if record.id==last_prot:
                  encontrada=True
            elif encontrada:
                beginning=time.time()
                result_handle = NCBIWWW.qblast("blastp", self.__db, record.format("fasta"), entrez_query='Homo sapiens [organism]')
                save_file.write(result_handle.read()+"\n")
                end=time.time()
                print("A proteína %s já foi submetida ao blast e demorou %s segundos. "%(record.id,end-beginning))
        
        
        
        