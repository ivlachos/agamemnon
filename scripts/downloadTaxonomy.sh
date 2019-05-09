#!/bin/bash


if [ ! -d ../Taxonomy/ ]
	then
		mkdir ../Taxonomy/
fi

echo "*** [Downloading NCBI Taxonomy]"

if [ ! -f  ../Taxonomy//dead_nucl.accession2taxid ]
	then
		wget -P ../Taxonomy/ ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz
		gunzip ../Taxonomy/dead_nucl.accession2taxid.gz
fi

if [ ! -f  ../Taxonomy/dead_wgs.accession2taxid ]
	then
		wget -P ../Taxonomy/ ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz
		gunzip ../Taxonomy/dead_wgs.accession2taxid.gz
fi

if [ ! -f  ../Taxonomy/nucl_gb.accession2taxid ]
	then
		wget -P ../Taxonomy/ ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
		gunzip ../Taxonomy/nucl_gb.accession2taxid.gz
fi

if [ ! -f  ../Taxonomy/nucl_gss.accession2taxid ]
	then
		wget -P ../Taxonomy/ ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz
		gunzip ../Taxonomy/nucl_gss.accession2taxid.gz
fi

if [ ! -f  ../Taxonomy/nucl_wgs.accession2taxid ]
	then
		wget -P ../Taxonomy/ ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
		gunzip ../Taxonomy/nucl_wgs.accession2taxid.gz
fi

if [ ! -f  ../Taxonomy/rankedlineage.dmp ]
	then
		wget -P ../Taxonomy/ ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
		tar -zxvf ../Taxonomy/new_taxdump.tar.gz rankedlineage.dmp
		mv rankedlineage.dmp ../Taxonomy/
		rm ../Taxonomy/new_taxdump.tar.gz
fi

