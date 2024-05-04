# `bmibdcaller_empirical`

This repository contains code used for validation of IBD caller benchmarking and
optimization findings in emprical data sets. 

The pipeline is part of a broader project `bmibdcaller`, which also includes the
following repositories:
- [bmibdcaller\_simulations](https://github.com/bguo068/posseleff_simulations) 
for simulations analysis.
- [ishare/ibdutils](https://github.com/bguo068/ishare)
for effienciently comparing two sets of inferred IBD segments.
- [tskibd](https://github.com/bguo068/tskibd) 
for obtaining true ibd segments from simulated genealogical trees. 


There are two main components in the repository, including
1. Scripts and notes to pre-process data available from [MalariaGEN Pf7](https://www.malariagen.net/resource/34/).
    - Data pre-filtering and downloading, inferring dominant alleles, and
    imputation. See notes [./input/Readme.md](./input/Readme.md).
    - Constructing different datasets. See notes [./datasets/Readme.md](./datasets/Readme.md).
2. The Nextflow pipeline to benchmark the performance of multiple IBD callers
before and after IBD caller parameter optimization: 
    - [Pipeline entry](./main.nf).

For the emprical analysis pipeline, the software environment and result folder
structure of this pipeline are similiar to the
[bmibdcaller\_simulations](https://github.com/bguo068/bmibdcaller_simulations)
repository. Details can be found Follow the [readme](https://github.com/bguo068/bmibdcaller_simulations/blob/main/README.md) from the simulation repository.


Examples of running the pipeline can be found in  
- [runs/r240113_opt/run.sh](runs/r240113_opt/run.sh).
- [runs/r240114_noopti/run.sh](runs/r240114_noopti/run.sh).