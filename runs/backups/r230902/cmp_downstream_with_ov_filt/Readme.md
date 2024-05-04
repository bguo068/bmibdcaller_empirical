# why
rerun this after re-downloading and re-filtering the samples with Fws>=0.95.
Sample size is smaller to than previous pipeline runs

# pipeline 

bmibdcaller_empirical commit:  c17166ffd26b1d12d7e6bc790f2e4005237f0861

# command
nextflow run /local/chib/toconnor_grp/bing/bmibdcaller_empirical/main.nf -profile sge -resume --filt_ibd_by_ov true
