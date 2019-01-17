# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 13:16:29 2019

@author: Joao
"""


from Bio import SeqIO

class FiltrarEssenciaisCodificantes:
    def __init__(self,genome_file,essential_file):
        self.__genome=SeqIO.read(genome_file,"genbank")
        self.__essential=essential_file
        self.__genes=[]
        self.__prot={}
        
    def getgenes(self):
        return self.__genes
    
    def print_genes(self):
        print(" ".join(self.__genes))
        
    def print_prot(self):
        for prot in self.__prot:
            print(">"+prot)
            print(self.__prot[prot]+"\n")
    
    def filtrar(self):
        genes_essenciais=open(self.__essential,"r")
        g_ess=[]
        lines=genes_essenciais.readlines()
        for gene in lines:
            g_ess.append(gene.strip("\n"))
        genes_essenciais.close()
        for i in range(len(self.__genome.features)):
            if self.__genome.features[i].type=="CDS":
                if "old_locus_tag" in self.__genome.features[i].qualifiers.keys():
                    num_tag=self.__genome.features[i].qualifiers["old_locus_tag"][0]
                    if num_tag in g_ess:
                        self.__genes.append(num_tag)
                        if "protein_id" in self.__genome.features[i].qualifiers.keys():
                            self.__prot[self.__genome.features[i].qualifiers["protein_id"][0]]=self.__genome.features[i].qualifiers["translation"][0]
                            
    def getprot(self):
        return self.__prot.items()
    
    def getseq(self,prot):
        return self.__prot[prot]
    
    def file_prot(self,file):
        new_file=open(file,"w")
        for prot in self.__prot:
            new_file.write(">"+prot+"\n")
            new_file.write(self.__prot[prot]+"\n \n")
        new_file.close()
    
    def verifica_gene(self,gene):
        if gene in self.__genes:
            return True
        else:
            return False
    
    def verifica_proteina(self,prot):
        if prot in self.__prot:
            return prot
        else:
            return None
     
"""       
featcds_essenciais=[]
g_ess=[]
genome=SeqIO.read("genesfull.gb","genbank")
genes_essenciais=open("genes_essenciais.txt","r")
lines=genes_essenciais.readlines()
new_genes_essenciais=open("genes_teste.fasta","w")
for gene in lines:
    g_ess.append(gene.strip("\n"))
genes_essenciais.close()
for i in range(len(genome.features)):
    if genome.features[i].type=="CDS":
        if "old_locus_tag" in genome.features[i].qualifiers.keys():
            num_tag=genome.features[i].qualifiers["old_locus_tag"][0]
            if num_tag in g_ess:
                featcds_essenciais.append(num_tag) #tags
                if "protein_id" in genome.features[i].qualifiers.keys():
                    new_genes_essenciais.write(">"+genome.features[i].qualifiers["protein_id"][0]+"\n")
                    new_genes_essenciais.write(genome.features[i].qualifiers["translation"][0]+"\n \n")
new_genes_essenciais.close()

    
print(len(featcds_essenciais))



from Bio.Blast import NCBIXML
from Bio import SeqIO
result_handle=open("resultado.xml")
blast_records=NCBIXML.parse(result_handle)
blasts=[]
for blast_record in blast_records: 
    blasts.append(blast_record.query)
blast_records.close()
print(blasts[len(blasts)-1])
print(len(blasts))

"""




