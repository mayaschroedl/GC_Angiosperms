#!/bin.sh
#From GFF3 file, we extract the needed information about exon and UTR to prepare the file that will be used on python to calcul GC content 
for SPECIES in $(cat list_species.txt)
	do
	echo ${SPECIES}
	sed '/#/d' ${SPECIES}.gene.gff3 > temp1
	sed 's/;/\t/g' temp1 > temp2
	cut -f1,3,4,5,7,8,9,10 temp2 > temp3 
	awk '{ print $7 "\t" $8 "\t" $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' temp3 > temp4
	sed '/gene/d'  temp4 > temp5
	sed '/mRNA[[:cntrl:]]/d' temp5 > temp6
	sed 's/ID=//g' temp6 > temp7
	sed 's/Parent=//g' temp7 > temp8
	cat head-exon temp8 > ${SPECIES}_exon_don.txt
	rm temp1 temp2 temp3 temp4 temp5 temp6 temp7 temp8
	done
exit