# sumstat2vcf

Sumstat2VCF is a wrapper script designed to efficiently convert Genome-Wide Association Study (GWAS) summary statistics into Variant Call Format (VCF) files. This tool simplifies the conversion process, making it faster and more manageable for downstream analyses. All necessary software components are dockerized, ensuring seamless execution without the need for additional installations. Furthermore, Sumstat2VCF allows for the seamless conversion of VCF files between different genome versions, such as hg19 and hg38.

## Prerequisites

* Docker: Ensure Docker is installed and running on your system.
* Input Files: Make sure you have the following files ready:
  * GWAS summary statistics file
  * Column dictionary file in JSON format
  * Reference genome FASTA file
  * dbSNP VCF file
  * Info file containing additional variant information (optional)
  * Target FASTA file (for liftover, optional)
  * Chain file (for liftover, optional)

## Running gwas2vcf

The command provided uses Docker to run sumstat2vcf. Below is the command to execute and an explanation of each parameter:

```
nohup docker run --rm \
    -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
    --user $(id -u):$(id -g) \
    jibinjv/gwas2vcf:v5 \
    python3 `pwd`/sumstat_to_vcf.py \
    --sumstat_file /mnt/disks/sdd/sumstat_folder/dbscan_clust_GenomicPCA_Correlation.N_weighted_GWAMA.results.txt.gz \
    --column_dict /mnt/disks/sdd/sumstat_folder/gwas_column.txt \
    --fasta /mnt/disks/sdd/sumstat_folder/human_g1k_v37.fasta \
    --dbsnp /mnt/disks/sdd/sumstat_folder/dbsnp.v153.b37.vcf.gz \
    --aliasfile NA \
    --gwas_outputname Cluster_GPCA_CORR \
    --output_folder /mnt/disks/sdd/sumstat_folder/vcf_files/ \
    --ncontrol 40000 \
    --ncase 30000 \
    --pvalue pvalue \
    --effectsize beta \
    --infofile /mnt/disks/sdd/sumstat_folder/variants.txt.gz \
    --infocolumn info \
    --eaffile /mnt/disks/sdd/sumstat_folder/variants.txt.gz \
    --eafcolumn af \
    --target_fasta /mnt/disks/sdd/sumstat_folder/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --chain_file /mnt/disks/sdd/sumstat_folder/GRCh37_to_GRCh38.chain.gz \
    --liftover Yes \
    > /mnt/disks/sdd/sumstat_folder/vcf_files/sumstat2vcf_nohup_log.txt 2>&1 & `
```

# Explanation of Command Parameters

* nohup: Runs the command in the background, allowing you to close the terminal without stopping the process.
* nohup: Runs the command in the background, allowing you to close the terminal without stopping the process.
* docker run --rm: Runs a Docker container and removes it after the command finishes.
* -v /mnt/disks/sdd/:/mnt/disks/sdd/: Mounts the host directory to the container directory.
* -v /home/jjohn41/:/home/jjohn41/: Mounts another host directory to the container directory.
* --user $(id -u):$(id -g): Runs the Docker container with the current user's ID and group ID to avoid permission issues.
* jibinjv/gwas2vcf:v5: Specifies the Docker image to use.
* python3 \pwd `/sumstat_to_vcf.py`: Path of the Python script.
* --sumstat_file: Specifies the GWAS summary statistics file with the full path.
* --column_dict: Specifies the column dictionary file in JSON format with the full path.
* --fasta: Specifies the reference genome FASTA file with the full path.
* --dbsnp: Specifies the dbSNP VCF file with the full path.
* --aliasfile: Alias file to map your GWAS summary stats chromosome name, default is "NA".
* --gwas_outputname: Specifies the name of the output GWAS file.
* --output_folder: Specifies the output folder with the full path.
* --ncontrol: Total number of controls in the study if binary or total sample size if continuous.
* --ncase: Total number of cases in the study if binary.
* --pvalue: Specifies whether the p-value is a regular p-value or -log10 p-value, choices are "mlogp" or "pvalue".
* --effectsize: Specifies whether the effect size is beta or odd ratio, choices are "beta" or "or".
* --infofile: Specifies the file with imputation info score, with the first four columns as chr, pos, a1, a2.
* --infocolumn: Specifies the info column name in the info file.
* --eaffile: Specifies the file with imputation eaf allele frequency, with the first four columns as chr, pos, a1, a2.
* --eafcolumn: Specifies the eaf allele frequency column name in the eaf file.
* --target_fasta: Specifies the target genome reference FASTA file for liftover.
* --chain_file: Specifies the chain file for liftover. The chain file can be downloaded from Ensembl.
* --liftover: Specifies whether liftover needs to be performed or not, choices are "No" or "Yes".
