#! /usr/bin/env python3

# software environment:
# mainly follow the instructure from MalariaGEN pf7 data user guide:
# https://malariagen.github.io/parasite-data/pf7/Data_access.html#variant-calls

# conda create -n pf7 python=3.10
# conda activate pf7
# pip install -q --no-warn-conflicts malariagen_data

import numpy as np
import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
import allel

# silence some warnings
dask.config.set(**{"array.slicing.split_large_chunks": False})
import malariagen_data

pf7 = malariagen_data.Pf7()
pf7_metadata = pf7.sample_metadata()
pf7_metadata.to_csv("pf7_meta.tsv", sep="\t", index=None)


# basic filters for variants and samples
extended_variant_dataset = pf7.variant_calls(extended=True)

## 1. remove samples not passing qc
sample_pass_qc = pf7_metadata["QC pass"].values
ds_sample_pass = extended_variant_dataset.isel(samples=sample_pass_qc)

filter_pass = ds_sample_pass["variant_filter_pass"].values
var_is_snp = ds_sample_pass["variant_is_snp"].values
is_biallelic = ds_sample_pass["variant_numalt"].values == 1

ds_bi = ds_sample_pass.isel(variants=(filter_pass & var_is_snp & is_biallelic))

ac0 = (ds_bi["call_genotype"].data[:, :, :2] == 0).sum(axis=[1, 2]).compute()
ac1 = (ds_bi["call_genotype"].data[:, :, :2] == 1).sum(axis=[1, 2]).compute()
nsamples = ds_bi.dims["samples"]
nploidy = ds_bi.dims["ploidy"]

non_single_double = (ac0 >= 2) & (ac1 >= 2)
smiss_le_0_8 = (1 - (ac0 + ac1) / (nsamples * nploidy)) <= 0.8

(non_single_double & smiss_le_0_8).mean()

# remove singleton and doubleton and sites with high missiningness
ds_bi = ds_bi.isel(variants=(non_single_double) & smiss_le_0_8)

ds_bi["call_AD"][:, :, :2]
# variant_AF/AC are zeros (not calculated)
# (biallelic_dataset["variant_AF"][:, 1].values > 0).sum()
# biallelic_dataset["variant_AC"][:, 0] >= 2
ds_bi["call_AD"].data[:, :, :2]
ds_bi["call_genotype"].data[:, :, :2]
ds_bi["variant_allele"].data[:, :2]


np.save("pos.npy", ds_bi["variant_position"].values)
np.save("chr.npy", ds_bi["variant_chrom"].values)
np.save("samples.npy", ds_bi["sample_id"].values)
np.save("alleles.npy", ds_bi["variant_allele"].values[:, :2])


# 	Array	Chunk
# Bytes	102.42 GiB	2.67 MiB
# Shape	(1696765, 16203, 2)	(7004, 100, 2)
# Data type	int16 numpy.ndarray
ad = ds_bi["call_AD"].data[:, :, :2].compute()
np.save("/data/bing/pf7_ad.npy", ad)


# # 54G
gt = ds_bi["call_genotype"].values
np.save("/data/bing/pf7_gt.npy", gt)
