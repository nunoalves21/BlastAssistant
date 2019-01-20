# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 17:17:39 2019

"""
from cmd import *
from FiltrarEssenciaisCodificantes import FiltrarEssenciaisCodificantes
from My_Blast import My_Blast
from ProcurarNaoHomologos import ProcurarNaoHomologos
import re
from CsvReaderLocus import CsvReaderLocus
import urllib.error
from AlignmentReader import AlignmentReader

class TargetingShell(Cmd):
    intro = """Principais comandos:
=================================================================================================                   
|    1. filtra_genes_essenciais - identificar genes essenciais e codificantes                   |
|    2. guardar - guardar as proteínas associadas a genes essenciais em formato fasta           |
|    3. verifica_gene - verifica se um gene indicado é um gene essencial                        |
|    4. ver_genes - mostra os genes essenciais e codificantes                                   |
|    5. ver_prot - mostra as proteínas associadas aos genes essenciais                          |
|    6. blast - fazer blast das proteínas guardadas contra o genoma humano                      |
|    7. blast_completar - se algum blast ficou a meio, execute este comando                     |
|    8. procurar_nao_homologos - identifica os genes não homólogos e regista num ficheiro       |
|    9. guardar_sumario - guardar informação relativa aos genes essenciais                      |
|    10. procurar_gene - imprime informação sobre um gene associado a uma determinada proteína  |
|    11. sumario_alinhamento - imprime informação sobre um alinhamento                          |
|    12. menu - quando necessitar de visualizar o menu escreva este comando                     |
|                           (Para mais comandos e ajuda escreva help)                           |       
=================================================================================================
"""
    prompt = "Targeting> "
    
    def do_filtra_genes_essenciais(self,arg):
        """Para identificar os genes essenciais no genoma do organismo pretendido é necessário que introduza um ficheiro com o genoma completo do mesmo em formato genbank e um ficheiro de texto com os genes essenciais com o seguinte formato:\n
            [old_locus_tag]
            [old_locus_tag]
            (...)
            ou o ficheiro csv retirado do site da OGEE
"""
        try:
            lst=arg.split(" ")
            if len(lst)==2:
                complete_genome=lst[0] #ficheiro de formato genbank com o genoma completo
                essenciais=lst[1] #ficheiro de formato txt ou csv com os genes essenciais
                rgx=re.compile("\d*.csv")
                global filtrar
                if rgx.search(essenciais)!=None: #se for um ficheiro com formato csv
                    reader=CsvReaderLocus(essenciais,"locus.txt") 
                    reader.read()
                    essenciais="locus.txt"
                    if filtrar!=None:
                        filtrar=None
                    filtrar=FiltrarEssenciaisCodificantes(complete_genome,essenciais)
                    filtrar.filtrar_locus_tag()
                    filtrar.print_genes()
                else:
                    if filtrar!=None:
                        filtrar=None
                    filtrar=FiltrarEssenciaisCodificantes(complete_genome,essenciais)
                    filtrar.filtrar_old()
                    filtrar.print_genes()
            else:
                print("Número de argumentos errados.")
        except:
            print("Ocorreu um erro a filtrar os genes")
                
    def do_guardar(self,arg):
        """Se já filtrou os genes essenciais e se quiser guardar as sequências de proteínas resultado da tradução dos genes essenciais, use o comando guardar usando como parâmetro o nome do ficheiro fasta em que pretende guardar as sequências de proteínas.""" 
        try:
            lst=arg.split(" ")
            if len(lst)==1:
                if filtrar==None:
                    print("Filtre primeiro os genes essenciais.")
                else:
                    newfile=lst[0]
                    filtrar.file_prot(newfile)
            else:
                print("Número de argumentos errados.")
        except:
            print("Erro a guardar")
    
    def do_guardar_sumario(self,arg):
        """Se já filtrou os genes essenciais e pretende guardar informação sobre estes genes num ficheiro use este comando usando como argumento um ficheiro."""
        try:
            lst=arg.split(" ")
            if len(lst)==1:
                if filtrar!=None:
                    filtrar.escrever_sumario(lst[0])
                else:
                    print("Filtre os genes essenciais primeiro, por favor.")
            else:
                print("Número de argumentos errado.")
        except:
            print("Erro a guardar o sumário")
                
    def do_sumario_gene(self,arg):
        """Comando que recebe como argumento o "locus_tag" ou o "old_locus_tag" de um gene e imprime as informações relativas a esse gene."""
        try:
            lst=arg.split(" ")
            if len(lst)==1:
                if filtrar!=None:
                    filtrar.ver_sumario(lst[0])
                else:
                    print("É necessário que filtre os genes essenciais primeiro.")
        except:
            print("Erro a dar sumario de gene")
            
    def do_verifica_gene(self,arg):
        """Verifica se um gene se encontra nos genes essenciais filtrados: recebe como parâmetro o "old_locus_tag" ou "locus_tag" do gene"""
        try:
            lst=arg.strip("\n").split(" ")
            if len(lst)==1:
                if filtrar.verifica_gene(lst[0])==True:
                    print("O gene é um gene essencial e codificante.")
                else:
                    print("O gene não é essencial e codificante.")
            else:
                print("Número de argumentos errados.") 
        except:
            print("Erro a verificar")
    
    def do_ver_genes(self,arg):
        """Comando que imprime os genes essenciais"""
        try:
            if filtrar!=None:
                filtrar.print_genes()
            else:
                print("Filtre primeiro os genes, por favor!")
        except:
            print("Erro a mostrar genes.")
    
    def do_ver_prot(self,arg):
        """Comando que imprime as proteínas associadas aos genes essenciais."""
        try:
            if filtrar!=None:
                filtrar.print_prot()
            else:
                print("Filtre primeiro os genes, por favor!")
        except:
            print("Erro a mostrar proteínas.")
    
    def do_verifica_proteina(self,arg):
        """Verifica se uma proteína se encontra entre as proteínas associadas aos genes essenciais filtrados: recebe como parâmetro o "protein_id" da proteína"""
        try:
            lst=arg.split(" ")
            if len(lst)==1:
                prot=filtrar.verifica_proteina(lst[0])
                if prot!=None:
                    print(prot+": \n"+filtrar.getseq(prot))
                else:
                    print("A proteína não está associada a nenhum gene essencial.")
            else:
                print("Número de argumentos errados.")
        except:
            print("Erro a verificar se a proteína está a associada a algum gene essencial")
    
    def do_sair(self, arg):
        "Sair do programa Targeting: sair"
        print('Obrigado por ter utilizado o Targeting, espero que tenha sido útil!')
        try:
            return True
        except:
            print("Erro a sair")
    
    def do_blast(self,arg):
        """Comando que faz um blast e coloca o seu resultado num ficheiro: recebe como argumentos um ficheiro fasta com as proteínas que vão ser submetidas ao blast, o nome do ficheiro de output com a extensão ".xml" e a base de dados."""
        try:
            lst=arg.strip("\n").split(" ")
            if len(lst)==3:
                blast=My_Blast(lst[0],lst[1],lst[2])
                blast.make_blast()
            else:
                print("Número de argumentos errados.")
        except urllib.error.URLError:
            print("""Erro de URL, se o blast se encontra a meio, por favor execute o comando "blast_completar".""")
        except:
            print("Erro a executar o blast.")
            
    def do_blast_completar(self,arg):
        """Comando que continua um blast que não se encontra terminado. Recebe como argumentos um ficheiro fasta com as proteínas que vão ser submetidas ao blast e o nome do ficheiro de output com a extensão ".xml" e a base de dados."""
        try:
            lst=arg.strip("\n").split(" ")
            if len(lst)==3:
                blast=My_Blast(lst[0],lst[1],lst[2])
                blast.make_blast_completar()
        except urllib.error.URLError:
            print("""Erro de URL, se o blast se encontra a meio, por favor execute o comando "blast_completar".""")
        except:
            print("Erro a fazer o blast.")
            
            
    def do_procurar_nao_homologos(self,arg):
        """Procura no resultado do blast os genes não homólogos. Recebe como argumentos um ficheiro xml resultado de um blast e o nome de um ficheiro de output"""
        try:
            lst=arg.split(" ")
            if len(lst)==2:
                nao_homologos=ProcurarNaoHomologos()
                nao_homologos.procura_naohomologos(lst[0],lst[1])
                nao_homologos.print_nao_homologos()
            else:
                print("Número de argumentos errado.")
        except:
            print("Erro a procurar proteínas/genes não homólogos.")

    def do_procurar_gene(self,arg):
        """Comando que recebe como argumento um id de uma proteína e imprime informação relativa ao gene associado."""
        try:
            lst=arg.split(" ")
            if len(lst)==1:
                if filtrar!=None:
                    filtrar.procurar_gene(lst[0])
                else:
                    print("Filtre os genes essenciais primeiro")
            else:
                print("Número de argumentos errado.")
        except:
            print("Erro a procurar genes.")

    def do_sumario_alinhamento(self,arg):
        """Comando que recebe um ficheiro em formato "fas" resultado de um alinhamento múltiplo e imprime informações sobre o consensus, zona mais conservada, etc.""" 
        try:
            lst=arg.split(" ")
            if len(lst)==1:
                align=AlignmentReader(lst[0])
                print("%s \n"%align.getAlignment())
                print("Consensus: \n%s \n"%align.getConsenso())
                print("A percentagem de blocos iguais em todas as sequências é de %f. \n"%align.getPerc())
                print("A zona mais conservada é a seguinte: %s \n"%align.getConservedDomain())
                print("A zona conservada inicia-se na posição: %d \n"%align.findConservedDomain())
            else:
                print("Número de argumentos errado.")
        except:
            print("Erro ao sumariar o alinhamento.")
    
    def do_menu(self,arg):
            print("""Principais comandos:
=================================================================================================                   
|    1. filtra_genes_essenciais - identificar genes essenciais e codificantes                   |
|    2. guardar - guardar as proteínas associadas a genes essenciais em formato fasta           |
|    3. verifica_gene - verifica se um gene indicado é um gene essencial                        |
|    4. ver_genes - mostra os genes essenciais e codificantes                                   |
|    5. ver_prot - mostra as proteínas associadas aos genes essenciais                          |
|    6. blast - fazer blast das proteínas guardadas contra o genoma humano                      |
|    7. blast_completar - se algum blast ficou a meio, execute este comando                     |
|    8. procurar_nao_homologos - identifica os genes não homólogos e regista num ficheiro       |
|    9. guardar_sumario - guardar informação relativa aos genes essenciais                      |
|    10. procurar_gene - imprime informação sobre um gene associado a uma determinada proteína  |
|    11. sumario_alinhamento - imprime informação sobre um alinhamento                          |
|    12. menu - quando necessitar de visualizar o menu escreva este comando                     |
|                           (Para mais comandos e ajuda escreva help)                           |       
=================================================================================================
""")
            

if __name__ == '__main__':
   filtrar=None 
   TargetingShell().cmdloop()
