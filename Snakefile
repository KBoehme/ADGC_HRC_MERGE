""" ADGC 2.0 pipeline. Inputs are collections of HRC imputed data and final
results are QC'ed dosage vcf files (1-22 chromosomes), and PLINK bed, bim, fam

NAMES = ['ACT1', 'ADC1', 'A DC3', 'ADC5', 'ADC7', 'BIOCARD', 'EAS', 'LOAD', 'MIRAGE', 'OHSU', 'ROSMAP1', 'TARC1', 'UKS',
         'UMVUTARC2', 'WASHU1', 'WHICAP', 'ACT2', 'ADC2', 'ADC4', 'ADC6', 'ADNI', 'CHAP2', 'GSK', 'MAYO', 'NBB',
         'RMAYO', 'ROSMAP2', 'TGEN2', 'UMVUMSSM', 'UPITT', 'WASHU2']

Run it kinda like so (For real without the -n):
cd /fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/logs
snakemake -np --verbose -j 999 --cluster-config ../scripts/cluster.json -s ../scripts/Snakefile --cluster \
"sbatch \
--ntasks {cluster.ntasks} \
--time {cluster.time} \
--mem {cluster.mem} \
--job-name {cluster.name} \
--mail-type {cluster.mail-type} \
--mail-user {cluster.mail-user} \
--parsable" \
--cluster-status ../scripts/slurm_status.py
"""

import os
from pathlib import Path
from python_scripts import utils

# This is the shared folder location, however due to lack of space cant be used.
# ROOT = "/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/"
ROOT = "/fslhome/fslcollab192/compute"
WORKFLOW_NAME = "combine_workflow"

ORIGINAL_DATA_SOURCE = os.path.join(ROOT, "ADGC_Data/ADGC_HRC")
PROCESS_DATA_SOURCE = os.path.join(ROOT, "ADGC_HRC_COMBINE", "process", WORKFLOW_NAME)
PRE_COMBINE_PREFIX = os.path.join(PROCESS_DATA_SOURCE, "pre_combine")
POST_COMBINE_PREFIX = os.path.join(PROCESS_DATA_SOURCE, "post_combine")

FINAL_DATA_SOURCE = os.path.join(ROOT, "ADGC_HRC_COMBINE", "final")
EXTRA_FOLDER = os.path.join(FINAL_DATA_SOURCE, "auxiliary")

RESOURCES = os.path.join("/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu", "resources")
SCRIPTS = os.path.join(ROOT, "scripts", "python_scripts")

LOGS = os.path.join(PROCESS_DATA_SOURCE, "logs")

MERGE_BASE_DATASET = "adc1"
# Explicity set rules to run locally. All others will run on cluster w/ defaults
localrules: create_sample_file_from_fam_per_study, keep_one_sample_file, set_fam_to_missing
DATASET_SKIP_LIST = []
PLINK_EXT = ["bed", "bim", "fam"]

DATASETS = [name for name in os.listdir(ORIGINAL_DATA_SOURCE) if os.path.isdir(os.path.join(ORIGINAL_DATA_SOURCE, name)) and name not in DATASET_SKIP_LIST]
print(DATASETS)
CHROMOSOME = list(range(1,23))

###### Various utility functions #######################
def get_sample_files(wildcards):
    sample_files = []
    for ds in DATASETS:
        if ds.lower() != MERGE_BASE_DATASET.lower():
            sample_files.append(os.path.join(PRE_COMBINE_PREFIX, ds.upper(), "{}.sample".format(ds.lower())))
    return sample_files

def get_gen_files(wildcards):
    gen_files = []
    for ds in DATASETS:
        if ds.lower() != MERGE_BASE_DATASET.lower():
            gen_files.append(os.path.join(PRE_COMBINE_PREFIX, ds, "chr{}_filtered.bgen".format(wildcards.chr)))
    return gen_files

def qctools_combine_params(wildcards):
    combined_command = []
    for ds in DATASETS:
        if ds.lower() != MERGE_BASE_DATASET.lower():
            combined_command.extend(["-g", os.path.join(PRE_COMBINE_PREFIX, ds, "chr{}_filtered.bgen".format(wildcards.chr)),
            "-s", os.path.join(PRE_COMBINE_PREFIX, ds, "{}.sample".format(ds.lower()))])
    return ' '.join(combined_command)

def get_all_combined_sample_files(wildcards):
    return [os.path.join(POST_COMBINE_PREFIX, f"chr{chr}.sample") for chr in CHROMOSOME]

#####################################################################

rule all:
    """
        Output of this workflow:
            Chr1-22 vcf.gz files
            Plink bed, bim, fam files
            Corrected fam flag (file indicating rule ran)
    """
    input:
        expand(os.path.join(FINAL_DATA_SOURCE, "vcf", "adgc_hrc_combined_chr{chr}.vcf.gz.tbi"), chr=CHROMOSOME),
        expand(os.path.join(FINAL_DATA_SOURCE, "plink", "adgc_hrc_merged_qced.{ext}"), ext=PLINK_EXT),
        # os.path.join(POST_COMBINE_PREFIX, "plink_prep", "updated_missing_final_fam.txt")


############### Plink production #####################
# Yay we are using plink now!
rule set_fam_to_missing:
    """ Fam files are missing some data thats in covars. Set to -9 sex and status
    """
    input:
        plink_fam=os.path.join(FINAL_DATA_SOURCE, "plink", "adgc_hrc_merged_qced.fam")
    output:
        fixed_fam=os.path.join(FINAL_DATA_SOURCE, "plink", "adgc_hrc_merged_qced.fam"),
        log_out=os.path.join(POST_COMBINE_PREFIX, "plink_prep", "updated_missing_final_fam.txt")
    params:
        temp_fam=os.path.join(FINAL_DATA_SOURCE, "plink", "adgc_hrc_merged_qced.fam.tmp")
    shell:
        "awk '{{print $1,$2,$3,$4,-9,-9}}' {input.plink_fam} > {params.temp_fam}; \
        sleep 2; \
        rm {input.plink_fam}; \
        sleep 3; \
        mv {params.temp_fam} {output.fixed_fam}; \
        sleep 3; \
        touch {output.log_out}"


rule convert_to_plink:
    """ Convert filtered bgen to plink format. Must use plink2 as we are using
    bgen format files."""
    input:
        combined_chr_bgen=os.path.join(POST_COMBINE_PREFIX, "plink_prep", "adgc_hrc_all_combined.bgen"),
        sample_file=os.path.join(POST_COMBINE_PREFIX, "combined.sample")
    output:
        expand(os.path.join(FINAL_DATA_SOURCE, "plink", "adgc_hrc_merged_qced.{ext}"), ext=PLINK_EXT)
    params:
        plink_base=os.path.join(FINAL_DATA_SOURCE, "plink", "adgc_hrc_merged_qced")
    shell:
        "plink2 --bgen {input.combined_chr_bgen} \
        --sample {input.sample_file} \
        --make-bed \
        --out {params.plink_base}"

rule cat_bgens:
    """Concatenate all the bgens into one humongous combined.bgen, to be converted to one humongous plink format. """
    input:
        gens=lambda wildcards: [os.path.join(POST_COMBINE_PREFIX, "maf_filtered", "bgen", f"chr{chr}_maf_filtered_rsid_updated.bgen") for chr in CHROMOSOME]
    output:
        concatenated_gen=os.path.join(POST_COMBINE_PREFIX, "plink_prep", "adgc_hrc_all_combined.bgen")
    shell:
        "cat-bgen -g {input.gens} -og {output.concatenated_gen}"

############### VCF production #####################

rule tabix_vcf:
    """ tabix vcf files"""
    input:
        dosage_vcf=os.path.join(FINAL_DATA_SOURCE, "vcf", "adgc_hrc_combined_chr{chr}.vcf.gz")
    output:
        output_tbi=os.path.join(FINAL_DATA_SOURCE, "vcf", "adgc_hrc_combined_chr{chr}.vcf.gz.tbi")
    shell:
        "tabix -p vcf {input.dosage_vcf}"

rule process_vcf:
    """ Fix the crappy vcf that comes out of qctools, strip ridiculous ",chr"
    from id, add AF, and contig header.
    """
    input:
        vcf=os.path.join(POST_COMBINE_PREFIX, "vcf", "chr{chr}_maf_filtered.vcf.gz"),
    output:
        out_vcf=os.path.join(FINAL_DATA_SOURCE, "vcf", "adgc_hrc_combined_chr{chr}.vcf.gz")
    shell:
        "export BCFTOOLS_PLUGINS=/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/programs/bcftools-1.8/plugins/; \
        utils_script.py fix-vcf {input.vcf} {wildcards.chr} | \
        bcftools +fill-tags -- -t AF | \
        utils_script.py dosage | \
        bgzip -c > {output.out_vcf}"

rule convert_to_vcf:
    """ Simply convert to vcf using qctools. Produces a vcf with GP and """
    input:
        combined_chr_gen=os.path.join(POST_COMBINE_PREFIX, "updated_rsid", "chr{chr}_combined_updated_rsid.gen.gz"),
        sample_file=os.path.join(POST_COMBINE_PREFIX, "combined.sample")
    output:
        vcf=os.path.join(POST_COMBINE_PREFIX, "vcf", "chr{chr}_maf_filtered.vcf.gz")
    log:
        os.path.join(LOGS, "qctools_convert_vcf", "chr{chr}_vcf_convert.log")
    params:
        threshold=0.5001 # Lowest threshold QCtools allows.
    shell:
        "/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/programs/gavinband-qctool-ba5eaa44a62f/build/release/qctool_v2.0.1 \
            -g {input.combined_chr_gen} \
            -s {input.sample_file} \
            -threshold {params.threshold} \
            -assume-chromosome {wildcards.chr} \
            -og {output.vcf} \
            -log {log}"

############### MAF filter Combined Dataset #####################
rule create_bgens:
    """ Just create bgens."""
    input:
        combined_chr_gen=os.path.join(POST_COMBINE_PREFIX, "updated_rsid", "chr{chr}_combined_updated_rsid.gen.gz")
    output:
        bgen=os.path.join(POST_COMBINE_PREFIX, "maf_filtered", "bgen", "chr{chr}_maf_filtered_rsid_updated.bgen")
    log:
        os.path.join(LOGS, "qctools_create_bgens", "chr{chr}.log")
    shell:
        "/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/programs/gavinband-qctool-ba5eaa44a62f/build/release/qctool_v2.0.1 \
            -g {input.combined_chr_gen} \
            -assume-chromosome {wildcards.chr} \
            -og {output.bgen} \
            -log {log}"

rule update_rsid_gen_file:
    """ Update rsid to real one if available, otherwise use fill chr:pos:A1:A2
    """
    input:
        gen=os.path.join(POST_COMBINE_PREFIX, "maf_filtered", "chr{chr}_maf_filtered.gen.gz"),
        dbsnp_ref=os.path.join(RESOURCES, "grch37_2017_ref", "chrom{chr}_ref.txt.gz")
    output:
        out_gen=os.path.join(POST_COMBINE_PREFIX, "updated_rsid", "chr{chr}_combined_updated_rsid.gen.gz")
    run:
        utils.rs_update_gen(input.gen, input.dbsnp_ref, output.out_gen)

rule filter_low_maf_variants:
    """ Remove low MAF SNPs. Use assume-chromosome to get the vcf chromosome
    column to properly populate"""
    input:
        combined_chr_bgen=os.path.join(POST_COMBINE_PREFIX, "chr{chr}_combined.bgen"),
        low_maf_snps=os.path.join(POST_COMBINE_PREFIX, "filters", "chr{chr}_maf_filter.txt"),
    output:
        filtered_vcf=os.path.join(POST_COMBINE_PREFIX, "maf_filtered", "chr{chr}_maf_filtered.gen.gz")
    log:
        os.path.join(LOGS, "qctools_merged_snp_stats", "chr{chr}_filter_maf_qctools.log")
    shell:
        "/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/programs/gavinband-qctool-ba5eaa44a62f/build/release/qctool_v2.0.1 \
            -g {input.combined_chr_bgen} \
            -assume-chromosome {wildcards.chr} \
            -og {output.filtered_vcf} \
            -excl-variants {input.low_maf_snps} \
            -log {log}"

rule write_low_maf_snps_per_chromosome:
    """ Using awk and the output of qctools -snp-stats, grab out all snps which
    are below the MAF specified in our param

    The columns are as follows, with header lines prefixed with # (like a vcf)

    maf_filter.txt format
    SNPID rsid chromosome position ref alt
    1:13380 1:13380 1 13380 C G

    snp-stats.txt file format
    22 22:16050435 NA 16050435 T C NA 1 1 401.823 0.177493 0.999558 0.000441525 0.000441525 C T 0.0278606 0.0278606 0 0 0 200.823 0.177493 2.08167e-17 0 201

    header of snp-stats.txt file
    alternate_ids	rsid	chromosome	position	alleleA	alleleB
    comment	HW_exact_p_value	HW_lrt_p_value	alleleA_count	alleleB_count
    alleleA_frequency	alleleB_frequency	minor_allele_frequency
    minor_allele	major_allele	info	impute_info	missing_proportion
    A	B	AA	AB	BB	NULL	total
    """
    input:
        os.path.join(POST_COMBINE_PREFIX, "filters", "chr{chr}snp-stats.txt")
    output:
        os.path.join(POST_COMBINE_PREFIX, "filters", "chr{chr}_maf_filter.txt")
    params:
        low_maf_threshold=".0001"
    shell:
        "echo SNPID rsid chromosome position ref alt > {output} && \
        grep -v -e '^#' -e '^alternate_id' {input} | awk '{{if ($14 < {params.low_maf_threshold}) print $2,$2,$1,$4,$5,$6}}' >> {output}"

rule generate_snp_statistics:
    """ Using Qctools, calculate per snp stats and output to a file. Then these
    can be used to gather a list of snps for exclusion such as super low MAF
    """
    input:
        combined_chr_gen=os.path.join(POST_COMBINE_PREFIX, "chr{chr}_combined.bgen")
    output:
        chr_snp_stats=os.path.join(POST_COMBINE_PREFIX, "filters", "chr{chr}snp-stats.txt")
    log:
        os.path.join(POST_COMBINE_PREFIX, "filters", "chr{chr}snp-stats_qctools.log")
    shell:
        "/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/programs/gavinband-qctool-ba5eaa44a62f/build/release/qctool_v2.0.1 \
            -g {input.combined_chr_gen} \
            -snp-stats \
            -osnp {output.chr_snp_stats} \
            -log {log}"

############### Combine Datasets #####################
rule keep_one_sample_file:
    """ qctools_combine_variants produces a sample file for each chromosome
    Rename chr1 to combined.sample and then remove all other sample files.
    An example of how to require many files but only reference one of them!
    """
    input:
        combined_chr_sample=get_all_combined_sample_files
    output:
        os.path.join(POST_COMBINE_PREFIX, "combined.sample")
    shell:
        "mv -n {input.combined_chr_sample[0]} {output}"

rule qctools_combine_variants:
    """ Merge each chrom from each study into a combined vcf.
    http://www.well.ox.ac.uk/~gav/qctool_v2/documentation/examples/combining.html
    i.e. All chr1 from all studies as input.
    """
    input:
        base_gen_file=os.path.join(PRE_COMBINE_PREFIX, MERGE_BASE_DATASET.upper(), "chr{chr}_filtered.bgen"),
        base_sample_file=os.path.join(PRE_COMBINE_PREFIX, MERGE_BASE_DATASET.upper(), MERGE_BASE_DATASET+".sample"),
        samples=get_sample_files,
        gens=get_gen_files
    output:
        combined_chr_gen=os.path.join(POST_COMBINE_PREFIX, "chr{chr}_combined.bgen"),
        combined_chr_sample=temp(os.path.join(POST_COMBINE_PREFIX, "chr{chr}.sample"))
    log:
        os.path.join(LOGS, "qctools_combine_variants", "chr{chr}_qctool.log")
    threads: 8
    params:
        combine_chrom_param=qctools_combine_params
    shell:
        "/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/programs/gavinband-qctool-ba5eaa44a62f/build/release/qctool_v2.0.1 \
            -g {input.base_gen_file} \
            -s {input.base_sample_file} \
            {params.combine_chrom_param} \
            -og {output.combined_chr_gen} \
            -os {output.combined_chr_sample} \
            -log {log} \
            -threads {threads}"

############### Pre combine Datasets #####################

rule filter_low_info_dup_variants:
    """ Remove low info snps. Takes low input filter snps and sample file.
    IMPORTANT: This uses assume chromosome when converting to bgen, otherwise
    these bgens will be messed when trying to convert to plink."""
    input:
        gen_file=os.path.join(ORIGINAL_DATA_SOURCE, "{ds}", "chr{chr}.gen.gz"),
        filter_variants=lambda wildcards: os.path.join(PRE_COMBINE_PREFIX, f"{wildcards.ds}", "filters", f"{wildcards.ds.lower()}_info_filter.txt")
    output:
        filtered_gen=os.path.join(PRE_COMBINE_PREFIX, "{ds}", "chr{chr}_filtered.bgen")
    log:
        os.path.join(LOGS, "filter_lowinfo", "{ds}", "chr{chr}_qctools.log")
    threads: 4
    shell:
        "/fslhome/fslcollab192/fsl_groups/fslg_KauweLab/compute/ADGC_2018_combined/alois.med.upenn.edu/programs/gavinband-qctool-ba5eaa44a62f/build/release/qctool_v2.0.1 \
            -g {input.gen_file} \
            -assume-chromosome {wildcards.chr} \
            -excl-variants {input.filter_variants} \
            -log {log} \
            -og {output.filtered_gen}"

rule write_dup_and_low_info_variants:
    """ This rule writes a single variant filter file to be used by qctools.
    It writes low info variants as well as duplicates.

    qctool -g example.bgen -og subsetted.bgen -excl-variants <filename>
    The specified file should contain a list of variants. Currently this must
    be a text file with six named columns; the first four must be SNPID, rsid,
    chromosome, position, followed by columns containing the first and second
    alleles. The -compare-variants-by option control how variants are matched
    to this file - see the page on sorting data for more information on this
    option.

    Using the gnomad exome data vcf, filter out the less common alts for any
    duplicate alelles."""
    input:
        info_file=os.path.join(ORIGINAL_DATA_SOURCE, "{ds}", "{ds_lower}_info.txt"),
        gnomad_ref_data=os.path.join(RESOURCES, "gnomad", "gnomad_exomes_af_snps_only.txt")
    output:
        filter_variant_file=os.path.join(PRE_COMBINE_PREFIX, "{ds}", "filters", "{ds_lower}_info_filter.txt"),
        dups_file=os.path.join(EXTRA_FOLDER, "duplicates", "{ds}_{ds_lower}_dups.txt")
    params:
        low_filter_threshold=".3"
    shell:
        "{SCRIPTS}/info_dup_filter.py {input.info_file} {input.gnomad_ref_data} {params.low_filter_threshold} {output.filter_variant_file} {output.dups_file}"

rule create_sample_file_from_fam_per_study:
    """ Creates the sample files from fam files"""
    input:
        os.path.join(ORIGINAL_DATA_SOURCE, "{ds}", "_{ds_lower}.fam")
    output:
        os.path.join(PRE_COMBINE_PREFIX, "{ds}", "{ds_lower}.sample")
    run:
        utils.sample_from_fam(input, output)

# rule write_low_info_snps_per_study:
#     """ Write low info snp to file.
#
#     qctool -g example.bgen -og subsetted.bgen -excl-variants <filename>
#     The specified file should contain a list of variants. Currently this must
#     be a text file with six named columns; the first four must be SNPID, rsid,
#     chromosome, position, followed by columns containing the first and second
#     alleles. The -compare-variants-by option control how variants are matched
#     to this file - see the page on sorting data for more information on this
#     option.
#     Example:
#
#     SNPID rsid chromosome position ref alt
#     """
#     input:
#         os.path.join(ORIGINAL_DATA_SOURCE, "{ds}/{ds_lower}_info.txt")
#     output:
#         os.path.join(PROCESS_DATA_SOURCE, "{ds}/filters/{ds_lower}_info_filter_variants.txt")
#     params:
#         low_filter_threshold=".03"
#     # log:
#     #     os.path.join(LOGS, "write_low_info_snps_per_study", "{ds}_{ds_lower}.log")
#     shell:
#         # This awk command simply prints the chr position which is non unique.
#         # "awk '{{if ($7 < {params.low_filter_threshold}) print $1}}' {input} > {output}"
#         # This awk command prints full variant info that should be used with the -excl-variants param for qctools.
#         "echo SNPID rsid chromosome position ref alt > {output} && \
#         awk 'FNR > 1 {{if ($7 < {params.low_filter_threshold}) split($1,var,\":\"); print $1,$1,var[1],var[2],$2,$3}}' {input} >> {output}"
