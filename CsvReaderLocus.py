# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 12:09:41 2019

"""

import csv

class CsvReaderLocus:
    def __init__(self,file,out):
        """Contrutor que recebe um ficheiro csv carregado da base de dados OGEE e um ficheiro de output com o formato txt"""
        self.file=file
        self.file_locus=out
        
    def read(self):
        with open(self.file) as csv_file:
            csv_reader=csv.reader(csv_file, delimiter=',')
            line_count=0
            out=open(self.file_locus,"w")
            for row in csv_reader:
                if line_count!=0:
                    if row[5]=="Essential":
                        out.write(row[0]+"\n")
                line_count+=1
            out.close()
            csv_file.close()
            
            
        
            

