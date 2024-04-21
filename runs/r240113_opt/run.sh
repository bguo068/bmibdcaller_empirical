#! /usr/bin/env bash
set -eEx -o pipefail
source ~/conda_devel.sh
conda activate bmibdcaller_simulations

cd /local/scratch/bing/bmibdcaller_empirical/runs/r240113_opt

nextflow run /local/chib/toconnor_grp/bing/bmibdcaller_empirical/main.nf -profile hq -resume \
    --filt_ibd_by_ov false  --ibdne_mincm 2.0 --ifm_mincm 2.0 --resdir opti_nofiltov_ne2_ifm2

nextflow run /local/chib/toconnor_grp/bing/bmibdcaller_empirical/main.nf -profile hq -resume \
    --filt_ibd_by_ov false  --ibdne_mincm 4.0 --ifm_mincm 4.0 --resdir opti_nofiltov_ne4_ifm4
