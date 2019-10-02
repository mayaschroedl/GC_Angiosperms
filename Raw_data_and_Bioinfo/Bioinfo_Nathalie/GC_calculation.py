from __future__ import print_function 

from Bio.Alphabet import single_letter_alphabet 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.SeqIO.Interfaces import SequentialSequenceWriter 
from Bio import SeqIO

from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

import string
import numpy as np
import pandas as pd
import csv
import os
import sys
import re


# list_species=["Athaliana","Atrichopoda","Egrandis","Macuminata","Spolyrhiza","Osativa"]
# list_species=["Acoerulea","Acomosus","Ahalleri","Ahypochondriacus","Alyrata","Bdistachyon","Boleraceacapitata","BrapaFPsc","Bstacei","Bstricta","Cclementina","Cgrandiflora","Cpapaya","Crubella","Csativus","Csinensis","Dcarota","Esalsugineum","Fvesca","Gmax","Graimondii","Kfedtschenkoi","Klaxiflora","Lusitatissimum","Mesculenta","Mguttatus","Mtruncatula","Othomaeum","Phallii","Ppersica","Ptrichocarpa","Pvulgaris","Rcommunis","Sbicolor","Sitalica","Slycopersicum","Spurpurea","Stuberosum","Sviridis","Tcacao","Tpratense","Vvinifera","Zmarina","Zmays"]
#
# list_species=["Mdomestica","Pvirgatum"]    #plus  long car plus de scaffold


#Functions to call to get the final data file (SPECIES_?_complet, with ?=RNA/exon/intron)
#prepare1(list_species)
#prepare2(list_species)
#merge_RNA(list_species)  
#merge_exon(list_species)
#merge_intron(list_species)


### 1 - Prepare data

## RNA Data 
#file1= SPECIES_RNA_don.txt , outfile= SPECIES_RNA_Len.txt
#The function 'longueur1' returns the lenght sequence for each RNA
def    longueur1(file1, outfile):
    file = open(file1, "r+")
    file2=open(outfile, "a")              #new file for sequence lenghts results
    file2.write("id \t len \n")
    head_line=file.readline()                     #head-line is not take in account in the loop
    for line in file:                             #loop to read the file      
        subch=line.split()                         #split the line in substring 
        len=int(subch[5])-int(subch[4])          #calcule sequence lenghts, subch[5]=end and subch[4]=start
        file2.write(subch[0])                  #id
        file2.write("\t")
        file2.write(str(len))                  #len
        file2.write("\n")
    file2.close()
    file.close()
         
 

## CDS-UTR Data
#file1=SPECIES_exon_don.txt, outfile=SPECIES_exon_len.txt
#The function 'longueur2' returns the lenght sequences for each cds and UTR sequence
def    longueur2(file1, outfile):                          
    file = open(file1, "r+")
    file2=open(outfile, "a")                            #new file for sequence lenghts results
    file2.write("id \t len \t num \n")
    headline=file.readline()                           #head-line is not take in account in the loop
    for line in file:                       #loop to read the file
        subch=line.split()                           
        len=int(subch[5])-int(subch[4])          #lenght
        id=subch[0]
        List=id.split(".")                 #split the substring id according to delimiter "." 
        file2.write(subch[0])                  #id
        file2.write("\t")
        file2.write(str(len))                 #len
        file2.write("\t")
        file2.write(str(List[-1]))              #n° CDS ou UTR: take the last suffix
        file2.write("\n")
    file2.close()
    file.close()
        
#file1=SPECIES_exon_don.txt, outfile=SPECIES_exon_nb.txt
#This function 'nbexons' returns the number of exons and intro per genes  
def nbexons(file1, outfile):                            
    file=open(file1, "r+")
    file2=open(outfile,"a")
    file2.write("keygene \t nbexon \t nbi \n")
    headline=file.readline()
    line1=file.readline()            #read the first ine to keep in memory the first keygene
    subch1=line1.split()
    L_keygene=[subch1[1]]
    if subch1[2]=='CDS':
        nbexon=1
    else:
        nbexon=0
    while 1:        #this king of loop have been chosen to be able to write the resutlt for the last keygene
        line=file.readline()
        subch=line.split()
        if line=="":                            #if we are at the end of the file. Careful on infinte loop if there is no empty line at the end of the file
            #print(nbexon)           
            nbi=nbexon-1
            file2.write(L_keygene[0])                   
            file2.write("\t")
            file2.write(str(nbexon))
            file2.write("\t")
            file2.write(str(nbi))
            break 
            
        else: 
            if subch[1]==L_keygene[0]:                        #keygene
                if subch[2]=='CDS':                #prend en compte que CDS, pas UTR
                    nbexon+=1
                    #print(subch[1])
                    #print(nbexon)
            else:
                #print(subch[1])
                nbi=nbexon-1
                file2.write(L[0])                   
                file2.write("\t")
                file2.write(str(nbexon))
                file2.write("\t")
                file2.write(str(nbi))
                file2.write("\n")
                L_keygene[0]=subch[1]          #replace the previous keygene par the current one
                nbexon=0
                if subch[2]=='CDS':
                    nbexon=1
    file.close()
    file2.close()


### 2- GC Calcul

#Basic function which returns the GC value of the sequence put in argument
def GC_calcul(sequece):
    return (GC(sequece))

## RNA Data
#fichier: fasta ESPECE.cds (sequences déjà séparées) et retourne un nouveau fichier avec id
#file2=SPECIES_cds.fa, outfile2=SPECIES_RNA_GCval.txt
#This function 'programme1' returns GC values of RNA sequence
def programme1(file2, outfile2):
    file=open(outfile2, "a")
    file.write("id \t GC_total \t GC1 \t GC2 \t GC3 \n")
    for seq_record in SeqIO.parse(file2, "fasta"):
        file.write(seq_record.name)                         #id du mRNA (correspond à son Name) court
        file.write("\t")
        file.write(str(GC_calcul(seq_record.seq)))                  #GC_total
        for k in range(0,3):                                  
            positionk=seq_record.seq[k::3]                       #extract each nucleotide at the k position
            file.write("\t")
            file.write (str(GC_calcul(positionk)))                #GCn°k
        file.write("\n")
    file.close()


## CDS-UTR Data
#Fct renvoyant la liste des numéros des chromosomes
#file1= SPECIES_exon
#This functions 'list_chr' returns in a list the names of the chromosomes + scaffolds
def list_chr(file1):                            
    file=open(file1, "r")
    headline=file.readline()                   #
    line1=file.readline()                      #read the first line to keep in memory the first chr
    subch1=line1.split()
    chr1=subch1[3]
    L=[chr1]
    i=0
    for line in file:
        subch=line.split()
        if subch[3]!=L[i]:
                    chr=subch[3]
                    L.append(chr)
                    i+=1
    #print(len(L))
    return (L)    

#argument var must be a string
#This function 'phase_seq' returns the position list to calculate first GC1, then GC2 and finally GC3, no matter is the phase of the sequence.
def phase_seq(var):                 #vérifier si le 1er nucleotide est le position 
    if var=="0": 
        Phase=[0,1,2]
    elif var=="1":
        Phase=[2,0,1]
    elif var=="2":
        Phase=[1,2,0]
    else: 
        Phase="NA"                   #return NA for UTR
    return(Phase)


#fichier1=SPECIES_exon_don.txt, fichier2=SPECIES_cds.fa 
#This function 'programme2' returns GC values of cds and UTR sequences
def programme2(fichier1, fichier2,outfile):                 
    file3=open(outfile, "a")                 #new file
    file3.write("id \t GC_total \t GC1 \t GC2 \t GC3 \n")
    L=list_chr(fichier1)                               #call fct to list the chromosome/scaffold name
    record_dict = SeqIO.index(fichier2, "fasta")          #index the lines of the fasta file
    for l in range(0, len(L)):
        chrom=L[l]
        sequence=record_dict[chrom]                         #extract the sequence corresponding to chr/scaffold l
        file1 = open(fichier1)                         #open  SPECIES_exon_don for each chr of the list L
        temp=open("tempo.txt","a")           #create temporary file only with line corresponding to the chr l
        for line in file1:
            subch1=line.split()                       
            if subch1[3]==L[l]:                    #if the line belong to the current chr/scaffold
                temp.write(line)
        file1.close()                             #close SPECIES_exon_don
        temp.close()        
        
        with open("tempo.txt") as temp_file:          #open to read it this time 
            for line in temp_file:
                subch=line.split()                         
                id=subch[0]                      
                file3.write(id)
                n=int(subch[4])                   #sequence start 
                m=int(subch[5])                   #sequence end
                my_seq=sequence.seq[n:m]          #extract the sequence portion from fasta file fichier2
                if subch[6]=="-":                  
                    my_seq=my_seq[::-1]           #reverse the sequence 
                file3.write("\t")
                file3.write(str(GC_calcul(my_seq)))
                Phases=phase_seq(subch[7])            #call function for sequence phase 
                j=0
                if Phases=="NA":                      #if UTR
                    for j in range (0,3):
                        file3.write("\t")
                        file3.write("NA")
                else:
                    while j< len(Phases):                       #Calcul first GC1, then GC2, then GC3
                        positionk=my_seq[Phases[j]::3]         
                        j+=1
                        file3.write("\t")
                        file3.write(str(GC_calcul(positionk)))
                file3.write("\n")
        os.remove("tempo.txt")
    file3.close() 


## #3- Obtain final data for RNA and CDS-UTR

#list_species= list of the species needed
#This funtion 'prepare1' runs the needed functions to get the different files for RNA 
def prepare1(list_species):
    in1="_RNA_don.txt"
    in2=".cds.fa" 
    out1= "_RNA_Len.txt"
    out2="_RNA_GCval.txt"
    for i in range(len(list_species)):
        ESPECE=list_species[i]
        print(ESPECE)
        file1=ESPECE + in1
        file2= ESPECE+ in2
        outfile1=ESPECE+ out1
        outfile2=ESPECE+ out2
        longueur1(file1, outfile1)      #--> SPECIES_RNA_Len.txt
        programme1(file2, outfile2)     #--> SPECIES_RNA_GCval.txt

#list_species= list of the species needed
#This funtion 'prepare2' runs the needed functions to get the different files for CDS and UTR
def prepare2(list_species):
    in1="_exon_don.txt"
    in2=".fa"
    out1="_exon_len.txt"
    out2="_exon_nb.txt"
    out3="_exon_GCval.txt"
    for i in range(len(list_species)):
        ESPECE=list_species[i]
        print(ESPECE)
        file1=ESPECE + in1
        file2= ESPECE+ in2
        outfile1=ESPECE+ out1
        outfile2=ESPECE+ out2
        outfile3=ESPECE+ out3
        
        longueur2(file1, outfile1)             #--> SPECIES_exon_len.txt
        nbexons(file1, outfile2)               #--> SPECIES_exon_nb.txt  
        programme2(file1, file2, outfile3)     #--> SPECIES_exon_GCval.txt


#This function merges all RNA file in the final file containing all information SPECIES_RNA_complet.txt
def merge_RNA(list_species):
    in1="_RNA_don.txt"
    out1= "_RNA_Len.txt"
    out2="_RNA_GCval.txt"
    out3="_exon_nb.txt"
    out4="_RNA_complet.txt"
    for i in range(len(list_species)):
        ESPECE=list_species[i]
        file1=ESPECE + in1
        outfile1=ESPECE+ out1
        outfile2=ESPECE+ out2
        outfile3=ESPECE+ out3
        outfile4=ESPECE+ out4
        
        d1= pd.read_table(file1,sep='\t',  names=['id','keygene','type','nChr','start','end'], low_memory=False)
        d2= pd.read_table(outfile1,sep='\t', names=['id','len'], low_memory=False)
        d3= pd.merge(d1,d2, on='id')
        d4= pd.read_table(outfile2,sep='\t', names=['id','GC_total','GC1','GC2','GC3'], low_memory=False)
        d5= pd.merge(d3,d4, on='id')                          
        d6= pd.read_table(outfile3, sep='\t', names=['keygene','nbexon', 'nbi'], low_memory=False)
        d7= pd.merge(d5, d6, on ='keygene')
        np.savetxt(outfile4, d7, fmt='%s', delimiter='\t', newline='\n',header='id\t keygene\t type\t nChr\t start\t end\t len\t GC_total\t GC1\t GC2\t GC3\t nbexon\t nbi', comments='' )


#This function merges all CDS-UTR file in the final file containing all information SPECIES_exon_complet.txt
def merge_exon(list_species):
    in1="_exon_don.txt"
    out1="_exon_len.txt"
    out2="_exon_nb.txt"
    out3="_exon_GCval.txt"
    out4="_exon_complet.txt"
    for i in range(len(list_species)):
        ESPECE=list_species[i]
        print(ESPECE)
        file1=ESPECE + in1
        outfile1=ESPECE+ out1
        outfile2=ESPECE+ out2
        outfile3=ESPECE+ out3
        outfile4=ESPECE+ out4
        
        f1= pd.read_table(file1,sep='\t',  names=['id','keygene','type','nChr','start','end','strand','phase'], low_memory=False)     #ou préciser dtype on import
        f2= pd.read_table(outfile1,sep='\t', names=['id','len','num'], low_memory=False)
        f3= pd.merge(f1,f2, on='id')
        f4= pd.read_table(outfile2 ,sep='\t',names=['keygene','nbexon','nbi'], low_memory=False )
        f5= pd.merge(f3,f4)                            
        f6= pd.read_table(outfile3, sep='\t', names=['id','GC_total','GC1','GC2','GC3'], dtype={'id': str} ,low_memory=False)
        f7=pd.merge(f5,f6, on='id')
        np.savetxt(outfile4, f7, fmt='%s', delimiter='\t', newline='\n', header='id\t keygene\t type\t nChr\t start\t end\t strand\t phase\t len\t num\t nbexon\t nbi\t GC_total\t GC1\t GC2\t GC3', comments='')
        


## #4 Introns Data
#file=SPECIES_exon_complet.txt, outfile=SPECIES_intron_coord.txt
#This function returns intron positions
def coord_intron(file,outfile):             #ESPECE_exon_complet, _intron_coord
    file1=open(file, "r+")
    file2=open(outfile,"a")                 #new file
    file2.write("id\tkeygene\ttype\tnChr\tnumi\tstart\tend")
    num=0                                  #counter for intron item
    headline=file1.readline()                      #read headline so not take in acocunt in the loop
    line1=file1.readline()                       #first line to extract id of the first gene
    subch1=line1.split()
    L=[subch1[1]]               #save keygene in memory
    if int(subch1[11])>0:             #for the first line, if nunmber of intron superior to 0
        if subch1[2]=='CDS':
            num=1
            file2. write("\n")
            file2.write(subch1[1])        #id    
            file2.write("_intron_")       #id    
            file2.write(str(num))         #id
            file2.write("\t")
            file2.write(subch1[1])        #keygene
            file2.write("\t")  
            file2.write("intron \t")      #type
            file2.write(subch1[3])        #nChr
            file2.write("\t")
            file2.write(str(num))        #numeral intron
            file2.write("\t")
            file2.write(subch1[5])      #end of the CDS = start of the next intron
        
    for line in file1:
        subch=line.split()
        if subch[1]==L[0]:                        #if same keygene
            if int(subch[11])>0:                  #if gene has at least 1 intron
                if subch[2]=='CDS':                #count only CDS, not UTR

                    if int(subch[9])==1:             #If it's the first CDS of a gene, subch[9]= num : numeral exon
                        num+=1
                        file2. write("\n")
                        file2.write(subch[1])        #id    
                        file2.write("_intron_")       #id    
                        file2.write(str(num))         #id
                        file2.write("\t")
                        file2.write(subch[1])      #keygene
                        file2.write("\t")
                        file2.write("intron")      #type
                        file2.write("\t")
                        file2.write(subch[3])       #nChr
                        file2.write("\t")
                        file2.write(str(num))       #intron item
                        file2.write("\t")
                        file2.write(subch[5])       #end exon > start for intron
                        
                    
                    elif int(subch[9])>num and int(subch[9])<int(subch[10]):    #if n° d'exon >counter et <number of exons
                        file2.write("\t")
                        file2.write(subch[4])          #start exon> end for intron
                        
                        num+=1
                        file2. write("\n")
                        file2.write(subch[1])        #id    
                        file2.write("_intron_")       #id    
                        file2.write(str(num))         #id
                        file2.write("\t")
                        file2.write(subch[1])     
                        file2.write("\t")
                        file2.write("\t")
                        file2.write(subch[3])      
                        file2.write("\t")
                        file2.write(str(num))
                        file2.write("\t")
                        file2.write(subch[5])       #end exon > début de intron
                    
                    elif int(subch[9])>num and int(subch[9])==int(subch[10]):           #last CDS of a gène
                        file2.write("\t")
                        file2.write(subch[4])          #start exon> end for intron
                    
                        
        else:                            #Gene change
            if int(subch[11])>0:         #if new gene has at least 1 intron     
                L[0]=subch[1]           #new gene saved
                num=0
                if subch[2]=='CDS':                    
                    num=1
                    file2. write("\n")
                    file2.write(subch[1])        #id
                    file2.write("_intron_")       #id
                    file2.write(str(num))         #id
                    file2.write("\t")
                    file2.write(subch[1])        #keygene
                    file2.write("\t")
                    file2.write("intron")
                    file2.write("\t")
                    file2.write(subch[3])       
                    file2.write("\t")
                    file2.write(str(num))
                    file2.write("\t")
                    file2.write(subch[5])        #start exon> end for intron
    file1.close()
    file2.close()


#fichier1=SPECIES_cds.fa.txt, fichier2=SPECIES_intron_coord.txt, outfile=SPECIES_intron_GCval.txt
#This function 'programme3' returns GC values intron sequences
def programme3(fichier1, fichier2,outfile):                 
    file3=open(outfile, "a")
    file3.write("id\tGC_total\n")
    L=list_chr(fichier2)                                 #call fct to list the chromosome/scaffold name
    record_dict = SeqIO.index(fichier1, "fasta")          #index the lines of the fasta file
    for l in range(0, len(L)):
        chrom=L[l]
        sequence=record_dict[chrom]              #extract the sequence corresponding to chr/scaffold l
        file2 = open(fichier2)                  #open  SPECIES_exon_don for each chr of the list L
        temp=open("tempo.txt","a")              #create temporary file only with line corresponding to the chr l
        for line in file2:
            subch=line.split()
            if subch[3]==L[l]:            #if the line belong to the current chr/scaffold
                temp.write(line)
        file2.close()                      
        temp.close()        
        
        with open("tempo.txt") as temp_file:           #open to read it this time 
            for line in temp_file:
                subch=line.split()                        
                file3.write(subch[0])              #id
                file3.write("\t")
                a=int(subch[5])                   #sequence start 
                b=int(subch[6])                   #sequence end
                n=min(a,b)                         #take the minimum, in case the gene is read in reverse
                m=max(a,b)
                my_seq=sequence.seq[n:m]                #extract the sequence portion from fasta file fichier1
                file3.write(str(GC_calcul(my_seq)))       #GC calcul
                file3.write("\n")
        
        os.remove("tempo.txt")
    file3.close() 
    

#This funtion 'merge_intron' runs the needed functions to get the different files for CDS and UTR and then merge them in SPECIES_intron_complet.txt
def merge_intron(list_species):
    in1="_exon_complet.txt"
    in2=".fa"
    out1="_intron_coord.txt"
    out2="_intron_GCval.txt"
    out3="_intron_complet.txt"
    for i in range(len(list_species)):
        ESPECE=list_species[i]
        print(ESPECE)
        file1=ESPECE + in1
        file2= ESPECE+ in2
        outfile1=ESPECE+ out1
        outfile2=ESPECE+ out2
        outfile3=ESPECE+ out3
        
        coord_intron(file1, outfile1)
        programme3(file2, outfile1, outfile2)
        
        f1= pd.read_table(outfile1,sep='\t',  names=['type', 'numi','keygene','nChr','start','end'], low_memory=False)
        f2= pd.read_table(outfile2,sep='\t', names=['type','numi','keygene','GC_total'],low_memory=False)
        f3= pd.merge(f1,f2, on=['id'])
        np.savetxt(outfile3, f3, fmt='%s', delimiter='\t', newline='\n')


