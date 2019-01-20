# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 13:16:29 2019

"""


from Bio import SeqIO

class FiltrarEssenciaisCodificantes:
    def __init__(self,genome_file,essential_file):
        """Construtor que recebe como parâmetros um ficheiro em formato genbank com o genoma completo do organismo patogénico, outro ficheiro de texto com os genes essenciais."""
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
    
    def filtrar_old(self):
        """Se o ficheiro com os genes essenciais for de formato txt, então terá os old_locus_tag desses genes. Esta função retorna todos os genes essenciais e codificantes"""
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
    
    def filtrar_locus_tag(self):
        """Módulo que filtra os genes essenciais e codificantes pelo locus_tag."""
        genes_essenciais=open(self.__essential,"r")
        g_ess=[]
        lines=genes_essenciais.readlines()
        for gene in lines:
            g_ess.append(gene.strip("\n"))
        genes_essenciais.close()
        for i in range(len(self.__genome.features)):
            if self.__genome.features[i].type=="CDS":
                if "locus_tag" in self.__genome.features[i].qualifiers.keys():
                    num_tag=self.__genome.features[i].qualifiers["locus_tag"][0]
                    if num_tag in g_ess:
                        self.__genes.append(num_tag)
                        if "protein_id" in self.__genome.features[i].qualifiers.keys():
                            self.__prot[self.__genome.features[i].qualifiers["protein_id"][0]]=self.__genome.features[i].qualifiers["translation"][0]
                            
    def getprot(self):
        return self.__prot.items()
    
    def getseq(self,prot):
        return self.__prot[prot]
    
    def file_prot(self,file):
        """Módulo que recebe como parâmetro um ficheiro vazio de formato fasta. O módulo escreve neste ficheiro as proteínas associadas aos genes essenciais.""" 
        new_file=open(file,"w")
        for prot in self.__prot:
            new_file.write(">"+prot+"\n")
            new_file.write(self.__prot[prot]+"\n \n")
        new_file.close()
    
    def verifica_gene(self,gene):
        """Módulo que recebe como parâmetro o "locus_tag" ou o "old_locus_tag" de um gene e verifica se este é essencial."""
        if gene in self.__genes:
            return True
        else:
            return False
    
    def verifica_proteina(self,prot):
        """Módulo que recebe como parâmetro o "protein_id" de uma proteína e verifica se esta está associada a algum gene essencial."""
        if prot in self.__prot:
            return prot
        else:
            return None
    
    def ver_sumario(self,gene):
        """Recebe como parâmetro o "locus tag" ou o "old_locus_tag" da um gene e imprime o nome do gene, a sua função, o produto da tradução, o id da proteína, a sua referência na base de dados e a tradução."""        
        handle=self.__genome
        for i in range(len(handle.features)):
            if "locus_tag" in handle.features[i].qualifiers.keys() and handle.features[i].type=="CDS":
                if handle.features[i].qualifiers["locus_tag"][0]==gene:
                    print(gene)
                    if "gene" in handle.features[i].qualifiers.keys():
                        print("O nome do gene é: %s."%handle.features[i].qualifiers["gene"][0])
                    if "note" in handle.features[i].qualifiers.keys():
                        print("Nota: %s."%handle.features[i].qualifiers["note"][0])
                    if "product" in handle.features[i].qualifiers.keys():
                        print("O produto da sua tradução é: %s."%handle.features[i].qualifiers["product"][0])
                    if "protein_id" in handle.features[i].qualifiers.keys():
                        print("O id da proteína é: %s."%handle.features[i].qualifiers["protein_id"][0])
                    if "db_xref" in handle.features[i].qualifiers.keys():
                        print("A referência na base de dados é: %s."%handle.features[i].qualifiers["db_xref"][0])
                    if "translation" in handle.features[i].qualifiers.keys():
                        print("A sequência da proteína traduzida é: %s"%handle.features[i].qualifiers["translation"][0])
            if "old_locus_tag" in handle.features[i].qualifiers.keys() and handle.features[i].type=="CDS":
                if handle.features[i].qualifiers["old_locus_tag"][0]==gene:
                    if "locus_tag" in handle.features[i].qualifiers.keys():
                        if "gene" in handle.features[i].qualifiers.keys():
                            print("O nome do gene é: %s."%handle.features[i].qualifiers["gene"][0])
                        if "note" in handle.features[i].qualifiers.keys():
                            print("Nota: %s."%handle.features[i].qualifiers["note"][0])
                        if "product" in handle.features[i].qualifiers.keys():
                            print("O produto da sua tradução é: %s."%handle.features[i].qualifiers["product"][0])
                        if "protein_id" in handle.features[i].qualifiers.keys():
                            print("O id da proteína é: %s."%handle.features[i].qualifiers["protein_id"][0])
                        if "db_xref" in handle.features[i].qualifiers.keys():
                            print("A referência na base de dados é: %s."%handle.features[i].qualifiers["db_xref"][0])
                        if "translation" in handle.features[i].qualifiers.keys():
                            print("A sequência da proteína traduzida é: %s"%handle.features[i].qualifiers["translation"][0])
                   
                        
    def escrever_sumario(self,file):
        """Módulo que recebe como parâmetro um ficheiro de texto e escreve informações sobre os genes essenciais"""
        genes_essenciais=open(self.__essential,"r")
        g_ess=[]
        lines=genes_essenciais.readlines()
        for gene in lines:
            g_ess.append(gene.strip("\n"))
        genes_essenciais.close()
        new_file=open(file,"w")
        handle=self.__genome
        for i in range(len(handle.features)):
            if "locus_tag" in handle.features[i].qualifiers.keys() and handle.features[i].type=="CDS":
                if handle.features[i].qualifiers["locus_tag"][0] in g_ess:
                    new_file.write("%s \n"%handle.features[i].qualifiers["locus_tag"][0])
                    if "gene" in handle.features[i].qualifiers.keys():
                        new_file.write("O nome do gene é: %s.\n"%handle.features[i].qualifiers["gene"][0])
                    if "note" in handle.features[i].qualifiers.keys():
                        new_file.write("Nota: %s. \n"%handle.features[i].qualifiers["note"][0])
                    if "product" in handle.features[i].qualifiers.keys():
                        new_file.write("O produto da sua tradução é: %s.\n"%handle.features[i].qualifiers["product"][0])
                    if "protein_id" in handle.features[i].qualifiers.keys():
                        new_file.write("O id da proteína é: %s.\n"%handle.features[i].qualifiers["protein_id"][0])
                    if "db_xref" in handle.features[i].qualifiers.keys():
                        new_file.write("A referência na base de dados é: %s.\n"%handle.features[i].qualifiers["db_xref"][0])
                    if "translation" in handle.features[i].qualifiers.keys():
                        new_file.write("A sequência da proteína traduzida é: %s \n \n"%handle.features[i].qualifiers["translation"][0])
            if "old_locus_tag" in handle.features[i].qualifiers.keys() and handle.features[i].type=="CDS":
                if handle.features[i].qualifiers["old_locus_tag"][0] in g_ess:
                    new_file.write("%s \n"%handle.features[i].qualifiers["old_locus_tag"][0])
                    if "gene" in handle.features[i].qualifiers.keys():
                        new_file.write("O nome do gene é: %s.\n"%handle.features[i].qualifiers["gene"][0])
                    if "note" in handle.features[i].qualifiers.keys():
                        new_file.write("Nota: %s. \n"%handle.features[i].qualifiers["note"][0])
                    if "product" in handle.features[i].qualifiers.keys():
                        new_file.write("O produto da sua tradução é: %s.\n"%handle.features[i].qualifiers["product"][0])
                    if "protein_id" in handle.features[i].qualifiers.keys():
                        new_file.write("O id da proteína é: %s.\n"%handle.features[i].qualifiers["protein_id"][0])
                    if "db_xref" in handle.features[i].qualifiers.keys():
                        new_file.write("A referência na base de dados é: %s.\n"%handle.features[i].qualifiers["db_xref"][0])
                    if "translation" in handle.features[i].qualifiers.keys():
                        new_file.write("A sequência da proteína traduzida é: %s \n \n"%handle.features[i].qualifiers["translation"][0])
        new_file.close()
        
    def procurar_gene(self,prot):
        """Módulo que recebe como parâmetro o id de uma proteína e imprime informação sobre o gene associado"""
        handle=self.__genome
        for i in range(len(handle.features)):
            if "protein_id" in handle.features[i].qualifiers.keys() and handle.features[i].type=="CDS":
                if prot==handle.features[i].qualifiers["protein_id"][0]:
                    self.ver_sumario(handle.features[i].qualifiers["locus_tag"][0])
