#!/bin/bash


if [ ! -d ../Taxonomy/ ]
	then
		mkdir ../Taxonomy/
fi

echo "*** [Downloading NCBI Taxonomy]"

while [ ! -f  ../Taxonomy/dead_nucl.accession2taxid ]; do
	wget -P ../Taxonomy/ https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz
	gunzip ../Taxonomy/dead_nucl.accession2taxid.gz
done

while [ ! -f  ../Taxonomy/dead_wgs.accession2taxid ]; do
	wget -P ../Taxonomy/ https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz
	gunzip ../Taxonomy/dead_wgs.accession2taxid.gz
done

while [ ! -f  ../Taxonomy/nucl_gb.accession2taxid ]; do
	wget -P ../Taxonomy/ https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
	gunzip ../Taxonomy/nucl_gb.accession2taxid.gz
done

while [ ! -f  ../Taxonomy/nucl_wgs.accession2taxid ]; do
	wget -P ../Taxonomy/ https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
	gunzip ../Taxonomy/nucl_wgs.accession2taxid.gz
done



while [ ! -f  ../Taxonomy/rankedlineage.dmp ]; do
	wget -P ../Taxonomy/ https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
	tar -zxvf ../Taxonomy/new_taxdump.tar.gz rankedlineage.dmp
	tar -zxvf ../Taxonomy/new_taxdump.tar.gz nodes.dmp
	tar -zxvf ../Taxonomy/new_taxdump.tar.gz names.dmp
	mv rankedlineage.dmp ../Taxonomy/
	mv nodes.dmp ../Taxonomy/
	mv names.dmp ../Taxonomy/
	rm ../Taxonomy/new_taxdump.tar.gz
done
