This directory contains a set of assembly files released as part of Phytozome release (these files have been renamed) used for GC content calculation and GC values files.

Full annoted genome of 52 angiosperm species were available on Phytozome database. Files have been renamed , here SPECIES match with a short form of the species name. 
Example:  Arabidopsis thaliana → Athaliana


Files in the phytozome_data subdirectory:
1. SPECIES_gene.gff3
GFF3 format representation of all mRNA sequences (UTR, CDS). Genomic coordinates are relative to the reference.
Sequence in column 1
2.  SPECIES.fa
Nucleotide FASTA format of the current genomic assembly, masked for repetitive sequence by RepeatMasker  (here softmasked sequence is in lower case).
3.  SPECIES_cds.fa
Nucleotide FASTA format file of all gene coding sequences, with alternative splice variants


Files in the script subdirectory:
1. list_species.txt
List of the 52 angiosperms species by their short form name, used by bash script to prepare data file for GC calculation.  
2. info_species.txt
File containing information about genome size, taxonomic group and full species names. 
3. prepa_RNA.sh
Bash script extracting information needed concerning RNA to create data file for GC calculs with python. Use list_species.  
4. prepa_exon.sh
Bash script extracting information needed concerning exon and UTR to create data file for GC calculs with python.  Use list_species.  
5. GC_calculation.py
Python script which calcules GC values for RNA, exons, UTR and introns and creates file containing all information obtained (SPECIES_?_complet).


Files in the GC_values subdirectory:
1. SPECIES_RNA_complet
Text file containing all information concerning RNA sequences.
	'id' matches to 'ID' from the GFF3 file for CDS and UTR sequences.
	'keygene' matches to 'Parent' from the GFF3 file for CDS and UTR sequences. 
	'nChr' is the name of the chromosome or scaffold on which is the sequence.
	'nbexon' is the number of exons per gene.
	'nbi' is the number of introns per gene.
2. SPECIES_exon_complet
Text file containing all information concerning CDS and UTR sequences.
	'id' matches to 'Name' from the GFF3 file for RNA sequences.
	'keygene' matches to 'ID' from the GFF3 file for RNA sequences. 
	'strand' defined as + (forward) or - (reverse) 
	'phase': One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on  
	'num': CDS item or UTR item
	'nbexon' is the number of exons per gene.
	'nbi' is the number of introns per gene.
3. SPECIES_intron_complet
	'id' of the intron has been invented by concatening the keygene + «_intron» + 'numi' value
	'keygene' of the intron matches to 'ID' of CDS from the GFF3 file for RNA sequences. 
	'numi' is the intron item
	
Extract from Read-me file of Phytozome:
Notes on JGI locus/gene naming convention: starting in 2013, JGI plant genome group began to use following naming pattern: 
   a) Prefix.NGn for stable chromosome scale genome assembly, example: Glyma.01G000100;
   b) Prefix.Zn for chromosome scale genome assembly, example: Eucgr.A00001, Pavir.Aa00001;
   c) Prefix.Nsn for fragmented genome assembly, example: Pahal.0001s0001
where N is chromosome number in a) or scaffold/contig number in c), Z is a letter representing a chromosome number like A for 1, B for 2 and so on in b), and n is locus/gene number on a chromosome or scaffold/contig. In both a) and b), a letter after last chromosome represents all scaffold/contig that are not mapped to a chromosome, for example, Glyma.U045300 (soybean has 20 chromosomes, 21st alphabet is U). Polyploid genome chromosome number can have a letter representing subgenome and hence N and Z could have one more letter like Pavir.Aa00001. Digit width of N in a) is variable ranging from 1 to 3, dot (period) is optional, and Prefix is made up from organism genus and species or chosen by community. Transcript name is locus name plus . plus digits (transcript number digit), for example, Glyma.01G000100.1. Initially transcript having digit 1 is longest but in subsequent gene annotation, transcript with digit 1 can be lost or not longest any more. The longest transcript should be looked up from GFF3 file with attribute longest=1. Initially, locus/gene number is ordered and increased by 100 from chromosome left to right in a), but this is very likely to change when genome assembly and/or gene annotation is updated. Please refer to external source for naming pattern for third party gene sets in Phytozome.

Nathalie Zeballos, 2e année Ingénieur Agronome, Montpellier Supagro, August 2018