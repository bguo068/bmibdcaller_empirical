# pipeline version
bmibdcaller_empirical: 8c9e825faf2db1986912f2d7486a2271d8728422
ibdutils: d9283c28d517ec2879315ea61c89581dea2f08bf


# Notes

- almost directly use scripts/pipeline from bmibdcaller_simulations pipeline
- no removal of IBD of TRMCA <1.5 generations

# run folders
```
/local/scratch/bing/bmibdcaller_empirical/runs/r230830
```

# command
```
nextflow run /local/chib/toconnor_grp/bing/bmibdcaller_empirical/main.nf -profile sge -resume
```
