# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 17:17:39 2019

@author: Joao
"""
from cmd import *
from FiltrarEssenciaisCodificantes import FiltrarEssenciaisCodificantes
from My_Blast import My_Blast
from ProcurarNaoHomologos import ProcurarNaoHomologos

class GingivShell(Cmd):
    intro = """Interpretador de comandos para o Targeting. Escrever help ou ? para listar os comandos disponíveis.\n
                                                Execute os seguintes comandos:
                =============================================================================================                    
                |    1. filtra_genes_essenciais - identificar genes essenciais                              |
                |    2. guardar - guardar as proteínas associadas a genes essenciais em formato fasta       |
                |    3. blast - fazer blast das proteínas guardadas contra o genoma humano                  |
                |    4. procurar_nao_homologos - identifica os genes não homólogos e regista num ficheiro   |
                |    5. verifica_gene - verifica se um determinado gene é essencial                         |
                ============================================================================================
                
                """
    prompt = 'Targeting> '
    
    def do_filtra_genes_essenciais(self,arg):
        """Para identificar os genes essenciais no genoma do organismo pretendido é necessário que introduza um ficheiro com o genoma completo do mesmo em formato genbank e um ficheiro de texto com os genes essenciais com o seguinte formato:\n
            [old_locus_tag]
            [old_locus_tag]
            (...) 
"""
        try:
            lst=arg.split(" ")
            if len(lst)==2:
                complete_genome=lst[0]
                essenciais=lst[1]
                global filtrar
                if filtrar!=None:
                    filtrar=None
                filtrar=FiltrarEssenciaisCodificantes(complete_genome,essenciais)
                filtrar.filtrar()
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
    
    def do_verifica_gene(self,arg):
        """Verifica se um gene se encontra nos genes essenciais filtrados: recebe como parâmetro o "old_locus_tag" do gene"""
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
        try:
            filtrar.print_genes()
        except:
            print("Erro a mostrar genes.")
    
    def do_ver_prot(self,arg):
        try:
            filtrar.print_prot()
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
        "Sair do programa Gingiv.: sair"
        print('Obrigado por ter utilizado o Gingiv, espero que tenha sido útil!')
        try:
            return True
        except:
            print("Erro a sair")
    
    def do_blast(self,arg):
        """Comando que faz um blast e coloca o seu resultado num ficheiro: recebe como argumentos um ficheiro fasta com as proteínas que vão ser submetidas ao blast e o nome do ficheiro de output com a extensão ".xml" e a base de dados."""
        try:
            lst=arg.strip("\n").split(" ")
            if len(lst)==3:
                blast=My_Blast(lst[0],lst[1],lst[2])
                blast.make_blast()
            else:
                print("Número de argumentos errados.")
        except:
            print("Erro a executar o blast.")
            
    def do_procurar_nao_homologos(self,arg):
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




if __name__ == '__main__':
   filtrar=None 
   GingivShell().cmdloop()
