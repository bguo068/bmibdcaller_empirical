#! /usr/bin/env bash
set -eEx -o pipefail

SCRATCH_DIR=/local/scratch/bing/bmibdcaller_empirical/runs/r240113_opt 
PIPELINE_DIR=/local/chib/toconnor_grp/bing/bmibdcaller_empirical

# source ~/conda_devel.sh
conda activate bmibdcaller_simulations

mkdir -p ${SCRATCH_DIR}
cd ${SCRATCH_DIR}

nextflow run ${PIPELINE_DIR}/main.nf -profile hq -resume \
    --filt_ibd_by_ov false  --ibdne_mincm 2.0 --ifm_mincm 2.0 --resdir opti_nofiltov_ne2_ifm2 && \
nextflow run ${PIPELINE_DIR}/main.nf -profile hq -resume \
    --filt_ibd_by_ov false  --ibdne_mincm 4.0 --ifm_mincm 4.0 --resdir opti_nofiltov_ne4_ifm4
