# AGAMEMNON
AGAMEMNON v0.1.0: an Accurate metagenomics and metatranscriptomics quantification analysis suite

### What is AGAMEMNON?
```AGAMEMNON is a comprehensive framework for the analysis of shotgun metagenomic/metatranscriptomic and host-specific RNA/DNA samples, enabling microbial quantification, differential abundance analyses and numerous exploratory visualizations.```

___

**NOTE:** Please download AGAMEMNON release v0.1.0 from the Releases tab!

___

### Dependencies
* Snakemake version >= 3.7
* R version >= 3.3
* Shiny R package >= 1.0.5
* Executable files for HISAT2 and Pufferfish are included in AGAMEMNON-v0.1.0 release, no installation required!
___

# Documentation
**BUILDING AN INDEX**</br>

The first step before start quantifying microbial abundances in NGS samples, is to build an index using the reference genomes of interest.</br>

To build an index, you need a multi-FASTA file containing microbial/viral genomes.</br>

**NOTE:** The headers of the FASTA file need to start with the NCBI Accession Number followed by anything else.</br>

For example:</br>
`>NC_013791.2 Bacillus pseudofirmus OF4, complete genome`</br>
`>NC_004061.1 Buchnera aphidicola str. Sg (Schizaphis graminum), complete genome`</br>
`>NC_002977.6 Methylococcus capsulatus str. Bath, complete genome`</br>

**The headers above, are all valid!**</br></br>
We also offer a set of ready-to-use microbial references (see papers Methods for further details)!</br>

To start building an index, navigate to the **../binaries/pufferfish/** directory and run the following command:</br>

**`./pufferfish index -r <reference.fasta> -o <output_directory> -p <num_of_threads>`**</br>

Were **<reference.fasta>** is a multi-FASTA file containing microbial/viral genomes, **<output_directory>** is the output directory and **<num_of_threads>** the number of threads to be allocated to Pufferfish.</br>

For additional indexing options (i.e k-mer size (default = 31)) simply run:</br>

**`./pufferfish index`**</br></br>

___


**QUANTIFYING MICROBIAL ABUNDANCES**</br>

In order to quantify the abundaces of microbial genomes using AGAMEMNON, you first need to fill the parameters of a configuration file (**config.yml**).</br>

The parameters of the **config.yml** file, are listed below:</br></br>

| Parameters | Description |
| :--- | :--- |
| BINARIES | Binary folder directory (default: ../agamemnon/binaries) - do not edit |
| CONFIG_DIR | A directory used internally by snakemake (default: ".") - do not edit |
| CONTROL_INDEX | Directory of the HISAT2 control index (i.e. spikeIn sequences, phix etc. - default: ../agamemnon/phix_index/phix) |
| HOST_INDEX | Directory of the HISAT2 host index, i.e. human, mouse etc. |
| INDEX_FASTQ | Directory of the index fastq file (for single-cell samples only) |
| PUFFERFISH_INDEX | Directory of the pufferfish index |
| RESULTS | Directory to write the results |
| SAMPLES_DIR | Directory of NGS samples |
| SCRIPTS | Directory of scripts (default: ../agamemnon/scripts) - do not edit |
| TAXONOMY_FILE | Directory of the NCBI taxonomy files (default: ../agamemnon/Taxonomy) - do not edit |
| ALGORITHM | Default: agamemnon - do not edit |
| CLEAR_ALL | Whether to clear everything except quantification results ("True") or not ("False") - default: "True" |
| FILES_EXT | Samples files extension, i.e. .fq, .fastq etc - (default: ".fq") |
| HOST_SAM | Whether to keep the SAM files produced by aligning reads with HISAT2 against host's genome/transcriptome ("True") or not ("False") - (default: "False") |
| MODE | Mode in which AGAMEMNON will be executed (1 or 2) - (default: 1)  |
| STRATEGY | "PE" for paired-end samples or "SE" for single-end samples - (default: PE) |
| TYPE | RNA or DNA host-specific samples, if DNA is selected, HISAT2 will be executed with the parameter --no-spliced-alignment (default: RNA) |
| MEM_MB | Default: 1 - (do not edit) |
| TPS | Number of threads to be used **per sample** - default: 8 |</br>

</br>**MODE 1:** Quantification of microbial abundances in **Shotgun Metagenomics/Metatranscriptomics** NGS samples.</br>

In MODE 1, AGAMEMNON will directly quantify the abundances of microbial genomes in Shotgun Metagenomics / Metatranscriptomics samples and thus, the parameters CONTROL_INDEX, HOST_INDEX, INDEX_FASTQ, HOST_SAM and TYPE are not necessary for the execution of the pipeline. Even if you chnage the default value ("NA"), it won't make any difference in the execution process.</br>

**MODE 2:** Quantification of microbial abundances in **host-specific RNA/DNA** NGS samples.</br>

<p align="justify">In MODE 2, at first, AGAMEMNON will align the sequencing reads against the host's genome/transcriptome using HISAT2, then against the selected control index (default: phix genome) and finally, it will quantify the abundances of microbial genomes in the remaining (unmapped) reads.</p>

To align sequencing reads against the host's genome and control sequences, AGAMEMNON is using HISAT2. You can download the ready-to-use indexes for both Homo sapiens (UCSC hg19) and Mus musculus (UCSC mm10) species following this link: http://www.ccb.jhu.edu/software/hisat/manual.shtml.</br>
Alternatively, you can build your own HISAT2 index using custom parameters and/or different genomes.</br>

**After filling the parameters above, the second step before the execution of AGAMEMNON is to navigate to the ../agamemnon/scripts/ directory and run the following command:**</br>

**`bash taxonomy.sh <reference.fasta>`**</br>

where the <reference.fasta> is the multi-FASTA file previously used to build the Pufferfish index. This command will both download the needed NCBI taxonomy files and find the TaxIDs for every Accession Number present in the multi-FASTA file.</br>

**Finally, to execute AGAMEMNON, use the following command:**</br>

**`snakemake --snakefile ../agamemnon/AGAMEMNON --cores <N> --resources mem_mb=<T>`**</br>

In the above command, "AGAMEMNON" is the snakemake file in the ../agamemnon/ directory, "--cores" refers to the maximum number of threads allowed to be used by AGAMEMNON and "mem_mb" refers to the maximum amount of RAM memory **(in MB)** allowed to be allocated to AGAMEMNON.</br>

Using the parameters "--cores", "mem_mb" and "TPS" (in the config.yml) file, AGAMEMNON will (automatically) make the best out of the available resources.</br>

If you want to display what would be done (test if the workflow is defined properly) but not actually execute AGAMEMNON, Snakemake gives the option **--dryrun** which prints a summary of the DAG of jobs (without really executing any job). This is useful in cases where a user has a large number of samples to analyse and she wants to be sure that everything is set up properly.</br>

After the execution of AGAMEMNON is finished, the results folder will contain the "../final/all" directory with a tab-delimited results file for every sample analysed, containing the microbial abundances together with their full lineage and TaxIDs. The results folder will also contain a folder for every sample, containing log files and if the value of "CLEAR_ALL" was setted to "False", it will also contain the **mappings.pam** file produced by Pufferfish and the **quant.sf** file produced by the quantification method. If AGAMEMNON was executed in MODE 2 with the "HOST_SAM" parameter setted to "True", the results folder will also contain the SAM files produced by HISAT2.</br>

___

**SHINY APPLICATION**</br>

We also offer a Shiny application where users can visualize and explore the quantification results produced by AGAMEMNON as well as perform differential abundance/expression analyses and conduct diversity index analysis.</br>

Before running the application, first you need to move the **/all folder** containing the quantification results to the **agamemnon/shinyApp** folder.</br>

In addition, you need to have a text file containing the phenotype information for your samples.</br>

The file must be named **phenotypes.tab** and needs to be a **tab-delimited** text file with the first column having the samples names and the rest of the columns having the phenotypic characteristics.</br>

The header of the first column must have the name **External ID**.</br>

**NOTE:** Samples names, are the file names inside the ../final/all directory without the .tab extension!

If you don't have phenotypic charasteristics for your samples, you need to provide at least a file including the first column with the samples names ("External ID") and a second column having a condition or random numbers. You can find an example of a phenotypes.tab file containing the minimum information inside the shinyApp folder.</br>

Having the **phenotypes.tab** file, you need to create a folder inside the /shinyApp directory with the name **phenoData** and move the phenotypes.tab file inside it.</br>

Now that you have the **../agamemnon/shinyApp/all/{sample_name}.tab** and **../agamemnon/shinyApp/phenoData/phenotypes.tab** files ready, open a terminal and run the following command:</br>

**`R -e "shiny::runApp('<full_path_to_shinyApp_folder>', launch.browser=TRUE)"`**</br>

This will open the Shiny Application using your default browser and you can start using the Application.</br>

You can also run the Shiny Application using R Studio. Just open the ../agamemnon/shinyApp/ui.R file with Rstudio and press the button **"Run App"** that appears in the top right corner.

