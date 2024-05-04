# pipeline version
bmibdcaller_empirical: ac46e356e18cbb138ef7d0eab98a7d2f3e4c153c
ibdutils: d9283c28d517ec2879315ea61c89581dea2f08bf


# Notes

- Use hmmIBD N_generation as TMRCA for segments called by hmmIBD
- Remove IBD from other called that overlapped with IBD segments from hmmIBD with TRMCA <1.5 generations

# run folders
```
/local/scratch/bing/bmibdcaller_empirical/runs/r230830
```

# command
```
nextflow run /local/chib/toconnor_grp/bing/bmibdcaller_empirical/main.nf -profile sge -resume --filt_ibd_by_ov true
```
