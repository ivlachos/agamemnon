# AGAMEMNON
AGAMEMNON version 0.1.0

### What is AGAMEMNON?
```AGAMEMNON is a comprehensive A-to-Z framework for the analysis of shotgun metagenomic/metatranscriptomic and host-specific RNA/DNA samples, enabling microbial quantification, differential abundance analyses and numerous visualizations.```

### Dependencies
* Snakemake version >= 3.7
* R version >= 3.3
* Shiny R package >= 1.0.5
* Executable files for Pufferfish, Cedar and HISAT2 are included in AGAMEMNON repository, no installation required.

### Documentation
The first step before start quantifying microbial genomes in NGS samples, is to build an index using the reference genomes of interest.</br></br>
To build an index, you need a multi-FASTA file containing microbial/viral genomes.</br><br>

**NOTE:** The headers of the FASTA file need to start with the genomes NCBI Accession Number followed by anything else.</br></br>

For example:</br>

`>NC_013791.2 Bacillus pseudofirmus OF4, complete genome`</br>
`>NC_004061.1 Buchnera aphidicola str. Sg (Schizaphis graminum), complete genome`</br>
`>NC_002977.6 Methylococcus capsulatus str. Bath, complete genome`</br></br>

**The headers above, are all valid!**</br></br>

