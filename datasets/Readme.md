# Goal

Based on the "clean" data prepared in the `input` folder, the goal is to prepare small data sets.
This includes two data sets for single-population (Eastern Southeast Asia and West
Africa), and one data set for a structured population (isolates from different
continents).

# Steps

1. Further clean the data by removing highly-related isolates
    - Based on IBD inferred via hap-ibd as it very fast for large sample size
    - The pipeline for inferring IBD can be found at
    [./01_call_hapibd/main.nf](./01_call_hapibd/main.nf)
    - Scripts to remove highly-related isoaltes based on IBD sharing:
    [./02_make_dataset/01_remove_related.py](./02_make_dataset/01_remove_related.py) for all samples and [](./02_make_dataset/03_remove_related_within_pop.py) for each population.
2. Construct data sets using the remaining isolates
    - Structured dataset 
        - Included populations: WSEA, ESEA, WAF, EAF and OCE.
        - Script: [./02_make_dataset/02_make_structured_subsets.py](./02_make_dataset/02_make_structured_subsets.py)
    - "Single population" datasets
        - Include a dataset for ESEA, and another for WAF.
        - Script: [./02_make_dataset/04_make_singlepop_subsets.py](./02_make_dataset/04_make_singlepop_subsets.py) 
    - Extract genotype data so that the resulting VCF files only contain targeted isolates
        - Script: [./02_make_dataset/05_subsetting_vcf.py](./02_make_dataset/05_subsetting_vcf.py)

Note: `./02_make_dataset/06_subsetting_vcf_reuse_posseleff_data.py` can be ignored, as it was intended
for exploratory analysis and was not included in the final analysis. 