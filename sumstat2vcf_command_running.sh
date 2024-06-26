
index=8
nohup docker run --rm \
    -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
    -v /home/jjohn41/:/home/jjohn41/ \
    --user $(id -u):$(id -g) \
    jibinjv/gwas2vcf:v5 \
    python3 `pwd`/sumstat_to_vcf.py \
    --sumstat_file /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/dbscan_clust_${index}_GenomicPCA_Correlation.N_weighted_GWAMA.results.txt.gz \
    --column_dict /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/gwas_column.txt \
    --fasta /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/human_g1k_v37.fasta \
    --dbsnp /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/dbsnp.v153.b37.vcf.gz \
    --aliasfile NA \
    --gwas_outputname Cluster_${index}_GPCA_CORR \
    --output_folder /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/vcf_files/ \
    --ncontrol NA \
    --ncase NA \
    --pvalue pvalue \
    --effectsize beta \
    --infofile /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/new_variants.txt.gz \
    --infocolumn info \
    --eaffile /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/new_variants.txt.gz \
    --eafcolumn af \
    --target_fasta /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --chain_file /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/GRCh37_to_GRCh38.chain.gz \
    --liftover Yes \
    > /mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/vcf_files/Cluster_${index}_GPCA_CORR_gwastovcf_nohup_log.txt 2>&1 &
