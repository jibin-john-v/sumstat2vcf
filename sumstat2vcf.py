import os,glob,argparse,json,re
import sys
import pandas as pd 
import numpy as np 
import concurrent.futures
import math


# Define command-line arguments
parser = argparse.ArgumentParser(description="This script will convert GWAS summary statistics to vcf format")

parser.add_argument('-sumstat_file', '--sumstat_file', help="Name of the sumstat file with full location", required=True)
parser.add_argument('-column_dict', '--column_dict', help="Name of the column in the sumstat file in JSON format", required=True)
parser.add_argument('-fasta', '--fasta', help="Genome reference FASTA file", required=True)
parser.add_argument('-dbsnp', '--dbsnp', help="dbSNP file", required=True)
parser.add_argument('-aliasfile', '--aliasfile', help=" alias file to map your GWAS summary stats chromosome name ",default="NA")
parser.add_argument('-gwas_outputname', '--gwas_outputname', help="Name of the GWAS file", required=True)
parser.add_argument('-output_folder', '--output_folder', help="Name of the output folder with full path", required=True)
parser.add_argument('-ncontrol', '--ncontrol', help="Total number of controls in the study if binary or total sample size if continuous", default="NA")
parser.add_argument('-ncase', '--ncase', help="Total number of cases in the study if binary", default="NA")
parser.add_argument('-pvalue', '--pvalue', help="Whether p value of -log10 p value", default="NA", choices=["mlogp", "pvalue"])
parser.add_argument('-effectsize', '--effectsize', help="Whether the effect size beta or odd ratio", default="beta", choices=["beta", "or"])
parser.add_argument('-infofile', '--infofile', help="file with imputation info score, first four columns should be chr,pos, a1,a2 ; as in the sumstat file", default="NA")
parser.add_argument('-infocolumn', '--infocolumn', help="Info column name in the infofile", default="infoscore")
parser.add_argument('-eaffile', '--eaffile', help="file with imputation eaf allele frequency, first four columns should be chr,pos, a1,a2 ; as in the sumstat file", default="NA")
parser.add_argument('-eafcolumn', '--eafcolumn', help="eaf allele frequency column name in the eaffile", default="eaf")
parser.add_argument('-target_fasta', '--target_fasta', help="Genome reference FASTA file of target for liftover", default="NA")
parser.add_argument('-chain_file', '--chain_file', help="chain file for liftover", default="NA")
parser.add_argument('-liftover', '--liftover', help="whether liftover need to to be perform or not", default="No",choices=["No", "Yes"])


args=parser.parse_args()

sumstat_file=args.sumstat_file
sumstat_column=args.column_dict
fasta=args.fasta
dbsnp=args.dbsnp
gwas_outputname=args.gwas_outputname
output_folder=args.output_folder
pvalue=args.pvalue
effectsize=args.effectsize
ncontrol=args.ncontrol
ncase=args.ncase
aliasfile=args.aliasfile
eaffile=args.eaffile
eafcolumn=args.eafcolumn
infofile=args.infofile
infocolumn=args.infocolumn
target_fasta=args.target_fasta
chain_file=args.chain_file
liftover=args.liftover


'''
sumstat_file="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/dbscan_clust_10_GenomicPCA_Correlation.N_weighted_GWAMA.results.txt.gz"
sumstat_column="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/gwas_column.txt"
fasta="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/human_g1k_v37.fasta" #Homo_sapiens_assembly38.fasta human_g1k_v37.fasta
dbsnp="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/dbsnp.v153.b37.vcf.gz" #dbsnp.v153.b37.vcf.gz dbsnp.v153.hg38.vcf.gz
gwas_outputname="IDP_Cluster10_Corr"
output_folder="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/vcf_files/"
pvalue="pvalue"
effectsize="beta"
infofile="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/new_variants.txt.gz"
infocolumn="info"
eaffile="/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/new_variants.txt.gz"
eafcolumn="af"
aliasfile="NA"
ncontrol="NA"
ncase="NA"
target_fasta='/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
chain_file='/mnt/disks/sdd/BrainImage_CLusterMetaSumstat/GPCA_Clusters/GRCh37_to_GRCh38.chain.gz'
liftover="Yes"
'''


n_parallel=8
sumstat_to_vcf_n_parallel=8
chunk_size=1300000

#Open the column json file
f = open(sumstat_column)
column_dict = json.load(f)


##constants 
required_columns=['chr_col','pos_col','snp_col','ea_col','oa_col','beta_col','se_col','pval_col','delimiter','header', 'build'] 
filt_column_dict= {key: value for key, value in column_dict.items() if value != 'NA'}
required_columns_in_sumstat=['chr_col','pos_col','snp_col','ea_col','oa_col','beta_col','se_col','pval_col']
required_columns_in_sumstat_dict={key: value for key, value in filt_column_dict.items() if key in required_columns_in_sumstat}
required_columns_in_json_dict={key: value for key, value in filt_column_dict.items() if key  in ['delimiter','header', 'build']}

# Check if all items in required_columns are present in filt_column_dict
if all(column in filt_column_dict for column in required_columns_in_sumstat):
    print("All items in required_columns are present in filt_column_dict.")
else:
    print("Not all items in required_columns are present in filt_column_dict.")
    sys.exit()  # Quit the program

#read sumstat file 
df=pd.read_csv(sumstat_file,sep=column_dict['delimiter'],engine="pyarrow")
df[column_dict['chr_col']] = df[column_dict['chr_col']].astype(str)
df[column_dict['chr_col']]=df[column_dict['chr_col']].replace('0X', 'X')
df[column_dict['pos_col']] = pd.to_numeric(df[column_dict['pos_col']], errors='coerce',downcast="integer")


##processing input sumstat
def process_sumstat_data(df, filt_column_dict, required_columns_in_sumstat_dict,ncontrol, ncase, effectsize, pvalue):
    # Process ncontrol_col
    if 'ncontrol_col' in filt_column_dict:
        required_columns_in_sumstat_dict['ncontrol_col'] = filt_column_dict['ncontrol_col']
        df[filt_column_dict['ncontrol_col']] = pd.to_numeric(df[filt_column_dict['ncontrol_col']], errors='coerce')
        df[filt_column_dict['ncontrol_col']] = df[filt_column_dict['ncontrol_col']].astype('int')
    elif 'ncontrol_col' not in filt_column_dict and ncontrol != "NA":
        df['ncontrol'] = ncontrol
        required_columns_in_sumstat_dict['ncontrol_col'] = 'ncontrol'
        filt_column_dict['ncontrol_col'] = 'ncontrol'
        df[filt_column_dict['ncontrol_col']] = pd.to_numeric(df[filt_column_dict['ncontrol_col']], errors='coerce')
        df[filt_column_dict['ncontrol_col']] = df[filt_column_dict['ncontrol_col']].astype('int')
    
    # Process 'ncase_col'
    if 'ncase_col' in filt_column_dict:
        required_columns_in_sumstat_dict['ncase_col'] = filt_column_dict['ncase_col']
        df[filt_column_dict['ncase_col']] = pd.to_numeric(df[filt_column_dict['ncase_col']], errors='coerce')
        df[filt_column_dict['ncase_col']] = df[filt_column_dict['ncase_col']].astype('int')
    elif 'ncase_col' not in filt_column_dict and ncase != "NA":
        df['ncase'] = ncase
        required_columns_in_sumstat_dict['ncase_col'] = 'ncase'
        filt_column_dict['ncase_col'] = 'ncase'
        df[filt_column_dict['ncase_col']] = pd.to_numeric(df[filt_column_dict['ncase_col']], errors='coerce')
        df[filt_column_dict['ncase_col']] = df[filt_column_dict['ncase_col']].astype('int')
    
    ##Process IFO and eaf
    rename_columns={'chr':column_dict['chr_col'],'pos':column_dict['pos_col'],'a1':column_dict['ea_col'],'a2':column_dict['oa_col']}
    merge_columns=[column_dict['chr_col'],column_dict['pos_col'],column_dict['ea_col'],column_dict['oa_col']]
    
    if eaffile != 'NA' or infofile != 'NA':
        if eaffile == infofile:
            # If eaffile and infofile are the same, read the file once
            info_df = pd.read_csv(eaffile, sep='\s+')
            info_df=info_df[['chr','pos','a1','a2',infocolumn, eafcolumn]]
            info_df['pos'] = pd.to_numeric(info_df['pos'], errors='coerce')
            info_df['chr'] = info_df['chr'].astype(str)
            info_df['chr']=info_df['chr'].replace('0X', 'X')
            info_df.rename(columns=rename_columns,inplace=True)
            df1 = pd.merge(df,info_df,
                left_on=[rename_columns['chr'], rename_columns['pos'], rename_columns['a1'], rename_columns['a2']],
                right_on=[rename_columns['chr'], rename_columns['pos'], rename_columns['a1'], rename_columns['a2']] )
            df2 = pd.merge(df,info_df,
                left_on=[rename_columns['chr'], rename_columns['pos'], rename_columns['a1'], rename_columns['a2']],
                right_on=[rename_columns['chr'], rename_columns['pos'], rename_columns['a2'], rename_columns['a1']] )
            df2=df2.drop(['A2_y','A1_y'],axis=1)
            df2=df2.rename(columns={"A1_x":"A1","A2_x":"A2"})
            df=pd.concat([df1,df2])
            df[[infocolumn, eafcolumn]] = df[[infocolumn, eafcolumn]].apply(lambda x: pd.to_numeric(x, errors='coerce'))
            required_columns_in_sumstat_dict['eaf_col'] = eafcolumn
            required_columns_in_sumstat_dict['imp_info_col'] = infocolumn
        else:
            # If eaffile and infofile are different, read them separately
            if eaffile != 'NA':
                eaf_df = pd.read_csv(eaffile, sep='\s+')
                eaf_df=eaf_df[['chr','pos','a1','a2',eafcolumn]]
                eaf_df['pos'] = pd.to_numeric(eaf_df['pos'], errors='coerce')
                eaf_df['chr'] = eaf_df['chr'].astype(str)
                eaf_df['chr']=eaf_df['chr'].replace('0X', 'X')
                eaf_df.rename(columns=rename_columns,inplace=True)
                df1 = pd.merge(df,eaf_df,
                    left_on=[rename_columns['chr'], rename_columns['pos'], rename_columns['a1'], rename_columns['a2']],
                    right_on=[rename_columns['chr'], rename_columns['pos'], rename_columns['a1'], rename_columns['a2']] )
                df2 = pd.merge(df,eaf_df,
                    left_on=[rename_columns['chr'], rename_columns['pos'], rename_columns['a1'], rename_columns['a2']],
                    right_on=[rename_columns['chr'], rename_columns['pos'], rename_columns['a2'], rename_columns['a1']] )
                df2=df2.drop(['A2_y','A1_y'],axis=1)
                df2[eafcolumn]=1-df2[eafcolumn]
                df2=df2.rename(columns={"A1_x":"A1","A2_x":"A2"})
                df=pd.concat([df1,df2])
                df[eafcolumn] = pd.to_numeric(df[eafcolumn], errors='coerce')
                required_columns_in_sumstat_dict['eaf_col'] = eafcolumn
            if infofile != 'NA':
                info_df = pd.read_csv(infofile, sep='\s+')
                info_df=info_df[['chr','pos','a1','a2',infocolumn]]
                info_df['pos'] = pd.to_numeric(info_df['pos'], errors='coerce')
                info_df['chr'] = info_df['chr'].astype(str)
                info_df['chr']=info_df['chr'].replace('0X', 'X')
                info_df.rename(columns=rename_columns,inplace=True)
                df=pd.merge(df,info_df,on=merge_columns,how="left")
                df[infocolumn] = pd.to_numeric(df[infocolumn], errors='coerce')
                required_columns_in_sumstat_dict['imp_info_col'] = infocolumn
    # Process 'imp_info_col'
    if infofile == 'NA'  and 'imp_info_col' in filt_column_dict:
        required_columns_in_sumstat_dict['imp_info_col'] = filt_column_dict['imp_info_col']
        df[filt_column_dict['imp_info_col']] = pd.to_numeric(df[filt_column_dict['imp_info_col']], errors='coerce')
    # Process 'eafcolumn'
    if eaffile == 'NA' and "eaf_col" in filt_column_dict:
        required_columns_in_sumstat_dict['eaf_col'] = filt_column_dict['eaf_col']
        df[filt_column_dict['eaf_col']] = pd.to_numeric(df[filt_column_dict['eaf_col']], errors='coerce')
    
    # Convert 'beta_col' to log scale if effectsize is "or"
    if effectsize == "or":
        df[filt_column_dict['beta_col']] = np.log(df[filt_column_dict['beta_col']])
    # Convert p-value to normal scale if pvalue is "mlogp"
    if pvalue == "mlogp":
        df[required_columns_in_sumstat_dict['pval_col']] = 10 ** (-df[required_columns_in_sumstat_dict['pval_col']])
    
    # Calculate 'imp_z_col'
    df['imp_z_col'] = df[filt_column_dict['beta_col']] / df[filt_column_dict['se_col']]
    required_columns_in_sumstat_dict['imp_z_col'] = 'imp_z_col'
    
    #correct the chromosome prefixe
    chr_column=filt_column_dict['chr_col']
    df[chr_column]=df[chr_column].astype('str').replace({"chr":""})
    df[chr_column]=df[chr_column].astype('str').replace({"23":"X","24":"Y","26":"MT"})
    return df, required_columns_in_sumstat_dict


df, required_columns = process_sumstat_data(df, filt_column_dict,required_columns_in_sumstat_dict, ncontrol, ncase, effectsize, pvalue)


#if the pvalue >=1, then convert to 0.99
df[required_columns_in_sumstat_dict['pval_col']]=np.where(df[required_columns_in_sumstat_dict['pval_col']]>=1,0.99,df[required_columns_in_sumstat_dict['pval_col']])


##Select the columns required for the analysis and convert data types to float 
def process_and_export_df(df, required_columns_in_sumstat_dict, output_folder, gwas_outputname, n_parallel, chunk_size=2000000):
    def save_to_tsv(chunk, output_file_name):
        chunk.to_csv(output_file_name, index=False,sep="\t")
    
    # Select required columns from DataFrame
    df2 = df[list(required_columns_in_sumstat_dict.values())]
    # Define float columns
    float_columns = [
        required_columns_in_sumstat_dict['beta_col'],
        required_columns_in_sumstat_dict['se_col'],
        required_columns_in_sumstat_dict['pval_col'],
        required_columns_in_sumstat_dict['imp_z_col'],
        required_columns_in_sumstat_dict['pos_col']
    ]
    # Convert selected columns to numeric
    df2[float_columns] = df2[float_columns].apply(pd.to_numeric, errors='coerce')
    num_chunks=math.ceil(df2.shape[0]/chunk_size)
    df3=np.array_split(df2,num_chunks)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as executor:
        futures = []
        for i in range(num_chunks):
            chunk = df3[i]
            output_file_name = f"{output_folder}{gwas_outputname}_part{i + 1}_vcf_input.tsv"
            future = executor.submit(save_to_tsv, chunk, output_file_name)
            futures.append(future)
        # Wait for all futures to complete
        for future in concurrent.futures.as_completed(futures):
            pass
    return df2, num_chunks

df2,num_chunks=process_and_export_df(df, required_columns_in_sumstat_dict, output_folder, gwas_outputname,n_parallel,chunk_size=chunk_size)




##Create json file which contain index number of columns 
def create_output_json(df2, required_columns_in_sumstat_dict, required_columns_in_json_dict, output_folder, gwas_outputname):
    df_columns = list(df2.columns)
    params = {}
    
    for key, value in required_columns_in_sumstat_dict.items():
        params[key] = df_columns.index(value)
    params = {**params, **required_columns_in_json_dict}
    with open(f"{output_folder}{gwas_outputname}.dict", "w") as outfile:
        json.dump(params, outfile)

create_output_json(df2, required_columns_in_sumstat_dict, required_columns_in_json_dict, output_folder, gwas_outputname)



##Run gwas2vcf
genome_build=required_columns_in_json_dict['build']
os.system(f'mkdir -p {output_folder}')

def parallel_sumstat_to_vcf(num_chunks, sumstat_to_vcf_n_parallel, output_folder, gwas_outputname, fasta, dbsnp, genome_build, aliasfile):
    def sumstat_to_vcf(sumstat_to_vcf_command):
        os.system(sumstat_to_vcf_command)
    
    def concatenate_and_process_vcf_files(output_folder, gwas_outputname, genome_build, sumstat_to_vcf_n_parallel):
        vcf_pattern = f"{output_folder}{gwas_outputname}_part*_build_{genome_build}_vcf.gz"
        vcf_file_list = sorted(glob.glob(vcf_pattern))
        err_pattern = f'{output_folder}{gwas_outputname}_gwastovcf_part*_errors.txt'
        err_file_list = sorted(glob.glob(err_pattern))
        
        # Define sorting key functions
        def sort_vcf_key(filename):
            match = re.search(r'part(\d+)_build', filename)
            return int(match.group(1)) if match else float('inf')
        
        def sort_err_key(filename):
            match = re.search(r'part(\d+)_errors', filename)
            return int(match.group(1)) if match else float('inf')
        
        # Sort the file lists using the sorting key functions
        vcf_files = sorted(vcf_file_list, key=sort_vcf_key)
        err_files = sorted(err_file_list, key=sort_err_key)
        
        # Join the sorted file names into strings
        vcf_files_str = ' '.join(vcf_files)
        err_files_str = ' '.join(err_files)
        
        # Concatenate the VCF files using bcftools
        output_vcf_file = f"{output_folder}{gwas_outputname}_build_{genome_build}_vcf.gz"
        os.system(f"bcftools concat --threads {sumstat_to_vcf_n_parallel} \
                    --allow-overlaps --rm-dups exact {vcf_files_str} | bgzip -c > {output_vcf_file}")
        os.system(f"tabix -p vcf {output_vcf_file}")
        
        # Remove the individual VCF files
        for vcf_file in vcf_files:
            os.remove(vcf_file)
        
        # Concatenate the error files
        output_err_file = f"{output_folder}{gwas_outputname}_gwastovcf_errors.txt"
        os.system(f"cat {err_files_str} > {output_err_file}")
        os.system(f"rm {err_files_str}")

    with concurrent.futures.ThreadPoolExecutor(sumstat_to_vcf_n_parallel) as executor:
        futures = []
        for i in range(num_chunks):
            input_file_name = f"{output_folder}{gwas_outputname}_part{i + 1}_vcf_input.tsv"
            output_file_name = f"{output_folder}{gwas_outputname}_part{i + 1}_build_{genome_build}_vcf"
            command = f'''python3 /app/main.py  \
                    --data {input_file_name} \
                    --ref {fasta} --dbsnp {dbsnp} \
                    --out {output_file_name} \
                    --id {gwas_outputname} \
                    --json {output_folder}{gwas_outputname}.dict '''
            
            if aliasfile != "NA":
                command += f'''  --alias {aliasfile} 2>{output_folder}{gwas_outputname}_gwastovcf_part{i + 1}_errors.txt'''
            else:
                command += f''' 2>{output_folder}{gwas_outputname}_gwastovcf_part{i + 1}_errors.txt'''
            
            future = executor.submit(sumstat_to_vcf, command)
            futures.append(future)
        
        # Wait for each job to complete and start the next one
        for future in concurrent.futures.as_completed(futures):
            pass
    
    concatenate_and_process_vcf_files(output_folder, gwas_outputname, genome_build, sumstat_to_vcf_n_parallel)
    
    # Remove intermediate TSV and TBI files
    sumstat_tsv_files = glob.glob(f"{output_folder}{gwas_outputname}_part*_vcf_input.tsv")
    tbi_files = glob.glob(f"{output_folder}{gwas_outputname}_part*_build_{genome_build}_vcf.gz.tbi")
    for file in sumstat_tsv_files + tbi_files:
        os.remove(file)




parallel_sumstat_to_vcf(num_chunks,sumstat_to_vcf_n_parallel,output_folder,gwas_outputname,fasta,dbsnp,genome_build,aliasfile)


def run_liftover_vcf(genome_build, output_folder, gwas_outputname, chain_file, target_fasta):
    if liftover == "Yes" and chain_file != "NA" and target_fasta != "NA":
        input_vcf = f"{gwas_outputname}_build_{genome_build}_vcf.gz"
        if genome_build == "GRCh38":
            # Convert to GRCh37
            output_vcf = f"{gwas_outputname}_build_GRCh37_vcf"
        elif genome_build == "GRCh37":
            # Convert to GRCh38
            output_vcf = f"{gwas_outputname}_build_GRCh38_vcf"
        # Run liftover_vcf.py in the background
        os.system(f'''python3 /usr/local/bin/liftover_vcf.py --fastafile {target_fasta} \
                      --chain_file {chain_file} --input_path {output_folder} --input_vcf {input_vcf} \
                      --output_vcf {output_vcf} --output_path {output_folder} >{output_folder}{output_vcf}_liftover_nohup_log.txt 2>&1 ''')


# Call the function with appropriate parameters
run_liftover_vcf(genome_build, output_folder, gwas_outputname, chain_file, target_fasta)