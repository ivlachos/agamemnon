# AGAMEMNON
AGAMEMNON version 0.1.0

### What is AGAMEMNON?
```AGAMEMNON is a comprehensive A-to-Z framework for the analysis of shotgun metagenomic/metatranscriptomic and host-specific RNA/DNA samples, enabling microbial quantification, differential abundance analyses and numerous visualizations.```

---

### Dependencies
* Snakemake version >= 3.7
* R version >= 3.3
* Shiny R package >= 1.0.5
* Executable files for Pufferfish, Cedar and HISAT2 are included in AGAMEMNON repository, no installation required.

___

### Documentation
The first step before start quantifying microbial genomes in NGS samples, is to build an index using the reference genomes of interest.</br>

To build an index, you need a multi-FASTA file containing microbial/viral genomes.</br>

**NOTE:** The headers of the FASTA file need to start with the genomes NCBI Accession Number followed by anything else.</br>

For example:</br>
`>NC_013791.2 Bacillus pseudofirmus OF4, complete genome`</br>
`>NC_004061.1 Buchnera aphidicola str. Sg (Schizaphis graminum), complete genome`</br>
`>NC_002977.6 Methylococcus capsulatus str. Bath, complete genome`</br>

**The headers above, are all valid!**</br></br>
We also offer a set of ready-to-use microbial references (see papers Methods for further information) which you can download following this link [link goes here]!</br>

**_BUILDING AN INDEX_**</br></br>
To start building an index, first navigate to the **../binaries/pufferfish/** directory and open the **config.json** file with a text editor.</br></br>

The values of the **config.json** file, are listed below:</br></br>

| Value | Description |
| :--- | :--- |
| ksize | The k-mer size to use (default: 27) |
| pufferfish | pufferfish directory (do not edit) |
| twopaco | twopaco directory (do not edit) |
| is_input_a_directory_to_fasta_files | whether you have a set of FASTA files (true) or a single multi-FASTA file (false) (default: false) |
| input_fasta | Directory to the multi-FASTA file |
| output_dir | Output directory to write the index |
| tmp_dir | Temporary folder to use while building the index (has to exist) |
| num_of_threads | Number of threads to use while building the index |
| twopaco_filter_size_microbiome | default: 37 (do not change) |
| twopaco_filter_size_humangenome | default: 36 (irrelevant for AGAMEMNON) |
| twopaco_filter_size_humantxome | default: 30 (irrelevant for AGAMEMNON) |
| twopaco_filter_size | default: 37 (do not change) |</br>

**Fill the values above, save the config.json file and start building an index using the following command:**</br></br>
**`bash index.sh config.json`**</br></br>
After the index building is done, you will end-up with 10 files in the output directory you previously set in the **config.json** file. Those files are: ctable.bin, eqtable.bin, info.json, mphf.bin, pos.bin, rank.bin, refAccumLenghts.bin, reflengths.bin, refseq.bin and seq.bin!</br>

**Please note** that the actual index directory will be a folder inside the output directory you set with a name following the pattern "fastafilename_ksize_fixed.puffidx"</br>

**_QUANTIFYING MICROBIAL GENOMES_**</br></br>

