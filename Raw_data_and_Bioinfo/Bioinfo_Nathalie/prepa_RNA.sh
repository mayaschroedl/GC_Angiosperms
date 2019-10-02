#!/bin.sh
#From GFF3 file, we extract the needed information RNA to prepare the file that will be used on python to calcul GC content 

for SPECIES in $(cat list_species.txt)
	do
	echo ${SPECIES}
	grep "mRNA[[:space:]]" ${SPECIES}.gene.gff3 > tempo1
	sed 's/;/\t/g' tempo1 > tempo2
	cut -f1,3,4,5,9,10 tempo2 > tempo3
	awk '{ print $6 "\t" $5 "\t" $2 "\t" $1 "\t" $3 "\t" $4}' tempo3 > tempo4
	sed 's/Name=//g' tempo4 > tempo5
	sed 's/ID=//g' tempo5 > tempo6
	cat head-RNA tempo6 > ${SPECIES}_RNA_don.txt
	rm tempo1 tempo2 tempo3 tempo4 tempo5 tempo6
	done
exit