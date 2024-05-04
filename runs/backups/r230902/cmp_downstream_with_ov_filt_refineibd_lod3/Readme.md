# Changes
1. set refined ibd lod to 3 (optized is 1.1 based on simulation)
2. add ESEA samples into the MP subworkflow
3. update plotting scripts

# pipeline version:

bmibdcaller_empirical commit: d6b8fe29bf8441f61cee5bb0d5d09df2e1af32ee

# command
nextflow run /local/chib/toconnor_grp/bing/bmibdcaller_empirical/main.nf -profile sge -resume --filt_ibd_by_ov true --refinedibd_lod 3
