This is the code for analysis of gain and loss of genes in whole-genome composition

It was intended for single-time computing, maybe for a few times only. Therefore, there is no interface and no guardrails.

First, the user needs to prepare the input file. The process starts with the selection of genomes. For instance, a dozen or two representative genomes from all or many branches of the animal kingdom. for each genome, get a file with as many gene models as possible. For most draft genomes gene models are predicted polypeptide sequences. All files should be in FASTA format. 

For each genome, make a two-letter mnemonic code (such as hs for Homo sapiens or am for Apis milifera, not necessarily an acronym, any two letters). Replace the name of each sequence in each FASTA file with the mnemonics, underscore, and sequence name, i.g. from ">genemodel123" to ">hs_genemodel123". 

Concatenate all FASTA files with gene models

Make BLAST indexed database

Run Blast search of all concatenated gene models against all concatenated gene models, with a tabulated output option. 

Now you have the input file for clustering. You can refine it by filtering out low e-value or low bitscore hits using select.c 

Use cluster.c to produce clusters of homologous sequences. This is a single-linkage clustering procedure linking all gene models connected through a common ancestor. 

Use map_sg.c to find out for each genome (represented by the initial set of gene models) which of the clusters are present and which are absent. The resulting file is a table where the rows are the clusters (loosely corresponding to gene families) and columns are the genomes. The last column lists a representative sequence. This file can be loaded into Excel and investigated with quick filters for the presence and absence of common ancestry gene clusters across genomes and taxa.

You can use fetch_samples.c to extract sequences of a particular cluster for further investigation.

