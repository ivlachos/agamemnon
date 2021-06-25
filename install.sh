#Pufferfish, Puffaligner and Cedar

echo ""
echo "Building Pufferfish, Puffaligner and Cedar ..."
echo ""

git clone https://github.com/COMBINE-lab/pufferfish.git
cd pufferfish
git checkout tags/v1.1.0
mkdir build
cd build
cmake ..
make


#HISAT2

echo ""
echo "Building HISAT2 ..."
echo ""

cd ../../
wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download
mv download hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip
rm hisat2-2.1.0-Linux_x86_64.zip


#Create the AGAMEMNON hierarchy and remove unnecessary files/folders

mkdir binaries
mv pufferfish binaries
cd binaries
mkdir cedar
cd pufferfish/build/src
cp pufferfish ../../
mv cedar ../../../cedar
cd ../../../../
mv hisat2-2.1.0 hisat2
mv hisat2 binaries


#Give executable permissions to python scripts
cd scripts
chmod +x accessions2taxids.py
chmod +x editYaml.py
chmod +x getTaxonomy.py
chmod +x indexSize.py
chmod +x single_cell-demultiplex.py

#END
