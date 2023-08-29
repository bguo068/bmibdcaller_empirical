# software environment
```
conda env create -f env.yml
```

# memory requirement
needs >250G for the step of retrieving the AD matrix from MalariaGEN

# scripts
```
01_download.py
02_infer_dmnt.py
03_filt_dom_gt.py
04_impute.py
```

# meta info
```
pf7_meta.tsv
```

# map by chromosome (use chr number instead of name)
```
impute/map_chr1.txt
impute/map_chr2.txt
impute/map_chr3.txt
impute/map_chr4.txt
impute/map_chr5.txt
impute/map_chr6.txt
impute/map_chr7.txt
impute/map_chr8.txt
impute/map_chr9.txt
impute/map_chr10.txt
impute/map_chr11.txt
impute/map_chr12.txt
impute/map_chr13.txt
impute/map_chr14.txt
```

# vcf by chromosome (use chr number, imputed, dominiant alleles)
```
impute/pf7_dom_0_9_imp_chr1.vcf.gz
impute/pf7_dom_0_9_imp_chr2.vcf.gz
impute/pf7_dom_0_9_imp_chr3.vcf.gz
impute/pf7_dom_0_9_imp_chr4.vcf.gz
impute/pf7_dom_0_9_imp_chr5.vcf.gz
impute/pf7_dom_0_9_imp_chr6.vcf.gz
impute/pf7_dom_0_9_imp_chr7.vcf.gz
impute/pf7_dom_0_9_imp_chr8.vcf.gz
impute/pf7_dom_0_9_imp_chr9.vcf.gz
impute/pf7_dom_0_9_imp_chr10.vcf.gz
impute/pf7_dom_0_9_imp_chr11.vcf.gz
impute/pf7_dom_0_9_imp_chr12.vcf.gz
impute/pf7_dom_0_9_imp_chr13.vcf.gz
impute/pf7_dom_0_9_imp_chr14.vcf.gz
```
