#! /usr/bin/env python3

# software environment:
# mainly follow the instructure from MalariaGEN pf7 data user guide:
# https://malariagen.github.io/parasite-data/pf7/Data_access.html#variant-calls

# conda create -n pf7 python=3.10
# conda activate pf7
# pip install -q --no-warn-conflicts malariagen_data=7.13.0 matplotlib

import numpy as np
import dask
import pandas as pd

# silence some warnings
dask.config.set(**{"array.slicing.split_large_chunks": False})
import malariagen_data

pf7 = malariagen_data.Pf7()
pf7_metadata = pf7.sample_metadata()
pf7_metadata.to_csv("pf7_meta.tsv", sep="\t", index=None)

# note this file is directly downloaded from the malariaGEN website instead of
# using the malariagen_data package
pf7_fws = pd.read_csv("./Pf7_fws.txt", sep="\t")

# merge with meta
metadata = pf7_metadata.merge(pf7_fws, on="Sample", how="left")


# basic filters for variants and samples
extended_variant_dataset = pf7.variant_calls(extended=True)

## 1. remove samples not passing qc
sample_pass_qc = metadata["QC pass"].values
is_monoclonal = metadata["Fws"] >= 0.95
sel_samples = sample_pass_qc & is_monoclonal
ds_sample_pass = extended_variant_dataset.isel(samples=sel_samples)

filter_pass = ds_sample_pass["variant_filter_pass"].values
var_is_snp = ds_sample_pass["variant_is_snp"].values
is_biallelic = ds_sample_pass["variant_numalt"].values == 1

ds_bi = ds_sample_pass.isel(variants=(filter_pass & var_is_snp & is_biallelic))

ac0 = (ds_bi["call_genotype"].data[:, :, :2] == 0).sum(axis=[1, 2]).compute()
ac1 = (ds_bi["call_genotype"].data[:, :, :2] == 1).sum(axis=[1, 2]).compute()
nsamples = ds_bi.dims["samples"]
nploidy = ds_bi.dims["ploidy"]

non_mac_le_10 = (ac0 >= 10) & (
    ac1 >= 10
)  # remove utra rare variants to reduce download size
smiss_le_0_8 = (1 - (ac0 + ac1) / (nsamples * nploidy)) <= 0.8

(non_mac_le_10 & smiss_le_0_8).sum()

# remove singleton and doubleton and sites with high missiningness
ds_bi = ds_bi.isel(variants=(non_mac_le_10) & smiss_le_0_8)

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
# Bytes	8.93 GiB	440.34 kiB
# Shape	(231596, 10348, 2)	(1281, 88, 2)
# Dask graph	77748 chunks in 5 graph layers
# Data type	int16 numpy.ndarray
ad = ds_bi["call_AD"].data[:, :, :2].compute()
np.save("/data/bing/pf7_ad.npy", ad)

# 	Array	Chunk
# Bytes	4.46 GiB	901.83 kiB
# Shape	(231596, 10348, 2)	(5247, 88, 2)
# Dask graph	12749 chunks in 4 graph layers
# Data type	int8 numpy.ndarray
gt = ds_bi["call_genotype"].values
np.save("/data/bing/pf7_gt.npy", gt)
