# pipeline version
bmibdcaller_empirical: e9a68f873eab78f436b9b32e95a8f38822e29150
ibdutils: d9283c28d517ec2879315ea61c89581dea2f08bf


# Notes

- Use hmmIBD N_generation as TMRCA for segments called by hmmIBD
- **Set TMRCA of hmmIBD IBD segments with length > 70% chrlen to 0** as N_generation < 1.5 missing a lot of tskibd IBD segment with TRMCA<1.5. This added step is an attempt to address the issue. See `misc/tmrca_hmmibd` folder for more details.
- Remove IBD from other called that overlapped with IBD segments from hmmIBD with TRMCA <1.5 generations

# run folders
```
/local/scratch/bing/bmibdcaller_empirical/runs/r230830
```

# command
```
nextflow run /local/chib/toconnor_grp/bing/bmibdcaller_empirical/main.nf -profile sge -resume --filt_ibd_by_ov true
```

# Key results

`runs/r230830/cmp_downstream_with_ov_filt/analysis_res/ne.png`
`runs/r230830/cmp_downstream_with_ov_filt2/analysis_res/ne.png`

This seems not changing too much of the results Ne (a dp around generation 20

