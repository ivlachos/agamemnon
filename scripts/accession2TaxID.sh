#!/bin/bash

REFERENCE_FASTA=$1
TAXONOMY_DIR=$2

mkdir $TAXONOMY_DIR/temp/

echo "[ Reading Accession Numbers from reference FASTA file... ]"
echo ""

grep ">" $REFERENCE_FASTA | awk '{print $1}' | awk -F ">" '{print $2}' > $TAXONOMY_DIR/accessions.txt
accessionsNum=$(wc -l $TAXONOMY_DIR/accessions.txt | awk '{print $1}')
nLines=0

echo "[ Finding Accession Numbers TaxIDs... ]"
echo ""

if [ $nLines -lt $accessionsNum ]
	then
		temp=0
		LC_ALL=C fgrep -w -f $TAXONOMY_DIR/accessions.txt $TAXONOMY_DIR/nucl_gb.accession2taxid | awk '{print $2"\t"$3}' > $TAXONOMY_DIR/temp/nucl_gb.tid
		temp=$(wc -l $TAXONOMY_DIR/temp/nucl_gb.tid | awk '{print $1}')
		nLines=$((nLines + temp))
fi

if [ $nLines -lt $accessionsNum ]
	then
		awk -F "\t" '{print $1}' $TAXONOMY_DIR/temp/nucl_gb.tid | grep -F -x -v -f - $TAXONOMY_DIR/accessions.txt > $TAXONOMY_DIR/restAccessions_gb.txt
		temp=0
		LC_ALL=C fgrep -w -f $TAXONOMY_DIR/restAccessions_gb.txt $TAXONOMY_DIR/dead_nucl.accession2taxid | awk '{print $2"\t"$3}' > $TAXONOMY_DIR/temp/dead_nucl.tid
		temp=$(wc -l $TAXONOMY_DIR/temp/dead_nucl.tid | awk '{print $1}')
		nLines=$((nLines + temp))
fi

if [ $nLines -lt $accessionsNum ] 
	then
		awk -F "\t" '{print $1}' $TAXONOMY_DIR/temp/dead_nucl.tid | grep -F -x -v -f - $TAXONOMY_DIR/restAccessions_gb.txt > $TAXONOMY_DIR/restAccessions_dead_nucl.txt
		temp=0
		LC_ALL=C fgrep -w -f $TAXONOMY_DIR/restAccessions_dead_nucl.txt $TAXONOMY_DIR/nucl_wgs.accession2taxid | awk '{print $2"\t"$3}' > $TAXONOMY_DIR/temp/nucl_wgs.tid
		temp=$(wc -l $TAXONOMY_DIR/temp/nucl_wgs.tid | awk '{print $1}')
		nLines=$((nLines + temp))
fi

if [ $nLines -lt $accessionsNum ]
	then
		awk -F "\t" '{print $1}' $TAXONOMY_DIR/temp/nucl_wgs.tid | grep -F -x -v -f - $TAXONOMY_DIR/restAccessions_dead_nucl.txt > $TAXONOMY_DIR/restAccessions_wgs.txt
		temp=0
		LC_ALL=C fgrep -w -f $TAXONOMY_DIR/restAccessions_wgs.txt $TAXONOMY_DIR/nucl_gss.accession2taxid | awk '{print $2"\t"$3}' > $TAXONOMY_DIR/temp/nucl_gss.tid
		temp=$(wc -l $TAXONOMY_DIR/temp/nucl_gss.tid | awk '{print $1}')
		nLines=$((nLines + temp))
fi

if [ $nLines -lt $accessionsNum ]
	then
		awk -F "\t" '{print $1}' $TAXONOMY_DIR/temp/nucl_gss.tid | grep -F -x -v -f - $TAXONOMY_DIR/restAccessions_wgs.txt > $TAXONOMY_DIR/restAccessions_gss.txt
		temp=0
		LC_ALL=C fgrep -w -f $TAXONOMY_DIR/restAccessions_gss.txt $TAXONOMY_DIR/dead_wgs.accession2taxid | awk '{print $2"\t"$3}' > $TAXONOMY_DIR/temp/dead_wgs.tid
		temp=$(wc -l $TAXONOMY_DIR/temp/dead_wgs.tid | awk '{print $1}')
		nLines=$((nLines + temp))
fi


#----------------------------------------------------------------
echo "[ Creating Accession Numbers to TaxIDs file... ]"
echo ""

cat $TAXONOMY_DIR/temp/*.tid > $TAXONOMY_DIR/accessionsTaxIDs.tab

echo "[ Removing unnecessary files... ]"
echo ""
rm $TAXONOMY_DIR/temp/*.tid


if [ -f $TAXONOMY_DIR/accessions.txt ]
	then
		rm $TAXONOMY_DIR/accessions.txt
fi

if [ -f $TAXONOMY_DIR/restAccessions_gb.txt ]
	then
		rm $TAXONOMY_DIR/restAccessions_gb.txt
fi

if [ -f $TAXONOMY_DIR/restAccessions_dead_nucl.txt ]
	then
		rm $TAXONOMY_DIR/restAccessions_dead_nucl.txt
fi

if [ -f $TAXONOMY_DIR/restAccessions_wgs.txt ]
	then
		rm $TAXONOMY_DIR/restAccessions_wgs.txt
fi

if [ -f $TAXONOMY_DIR/restAccessions_gss.txt ]
	then
		rm $TAXONOMY_DIR/restAccessions_gss.txt
fi

rm -r $TAXONOMY_DIR/temp/

