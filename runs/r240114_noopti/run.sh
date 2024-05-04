#! /usr/bin/env bash
set -eEx -o pipefail

SCRATCH_DIR=/local/scratch/bing/bmibdcaller_empirical/runs/r240114_noopti 
PIPELINE_DIR=/local/chib/toconnor_grp/bing/bmibdcaller_empirical

# source ~/conda_devel.sh
conda activate bmibdcaller_simulations

mkdir -p ${SCRATCH_DIR}
cd ${SCRATCH_DIR} 

nextflow run ${PIPELINE_DIR}/main.nf -profile hq -resume \
    --hapibd_minoutput 2.0 \
    --hapibd_minseed 2.0 \
    --hapibd_minextend 1.0 \
    --hapibd_maxgap 1000 \
    --hapibd_minmarkers 100 \
    --hmmibd_n 999999 \
    --hmmibd_m 5 \
    --isorelate_imiss  0.1 \
    --isorelate_vmiss  0.1 \
    --isorelate_min_snp  20 \
    --isorelate_minmaf  0.01 \
    --refinedibd_length 2.0 \
    --refinedibd_lod  3.0 \
    --refinedibd_scale  0 \
    --refinedibd_window  40.0 \
    --tpbwt_template_opts 0 \
    --tpbwt_Lm 300 \
    --tpbwt_Lf 2.0 \
    --tpbwt_use_phase_correction 1 \
    --filt_ibd_by_ov false  --ibdne_mincm 2.0 --ifm_mincm 2.0 --resdir noopti_nofiltov_ne2_ifm2 \
    && \
nextflow run ${PIPELINE_DIR}/main.nf -profile hq -resume \
    --hapibd_minoutput 2.0 \
    --hapibd_minseed 2.0 \
    --hapibd_minextend 1.0 \
    --hapibd_maxgap 1000 \
    --hapibd_minmarkers 100 \
    --hmmibd_n 999999 \
    --hmmibd_m 5 \
    --isorelate_imiss  0.1 \
    --isorelate_vmiss  0.1 \
    --isorelate_min_snp  20 \
    --isorelate_minmaf  0.01 \
    --refinedibd_length 2.0 \
    --refinedibd_lod  3.0 \
    --refinedibd_scale  0 \
    --refinedibd_window  40.0 \
    --tpbwt_template_opts 0 \
    --tpbwt_Lm 300 \
    --tpbwt_Lf 2.0 \
    --tpbwt_use_phase_correction 1 \
    --filt_ibd_by_ov false  --ibdne_mincm 4.0 --ifm_mincm 4.0 --resdir noopti_nofiltov_ne4_ifm4
