import numpy as np
import matplotlib.pyplot as plt
from subprocess import run

pos = np.load("pos.npy")
chrom = np.load("chr.npy", allow_pickle=True)
samples = np.load("samples.npy", allow_pickle=True)
alleles = np.load("alleles.npy", allow_pickle=True)
dgt_orig = np.load("/data/bing/pf7_gt_dominant_0_9.npy")
dgt = dgt_orig.copy()

# 1. remove sites of rare variants or high missiningness
nsites, nsamples = dgt.shape

af = (dgt == 1).sum(axis=1) / ((dgt == 0).sum(axis=1) + (dgt == 1).sum(axis=1))
no_rare_vars = (af > 0.01) & (1 - af >= 0.01)

site_miss = (dgt == -1).sum(axis=1) / nsamples
no_extreme_per_site_miss = site_miss <= 0.8

sel_sites = no_extreme_per_site_miss & no_rare_vars

pos = pos[sel_sites].copy()
chrom = chrom[sel_sites].copy()
alleles = alleles[sel_sites].copy()
dgt = dgt[sel_sites, :].copy()  # after filtering (25681, 16203)


# 2. remove samples with high missingness
nsites, nsamples = dgt.shape
sample_miss = (dgt == -1).sum(axis=0) / nsites
plt.hist(sample_miss)
# max missingness pers sample is 0.5 , no need to filter this round
sample_miss.max()

# 3. keep sites with <- 10% missing (2nd round more stringent)
site_miss = (dgt == -1).sum(axis=1) / nsamples
plt.hist(site_miss)
sel_sites = site_miss <= 0.1

pos = pos[sel_sites].copy()
chrom = chrom[sel_sites].copy()
alleles = alleles[sel_sites].copy()
dgt = dgt[sel_sites, :].copy()

# 4. keep samples with < 10% missingness (2nd round more stringent)
nsites, nsamples = dgt.shape
sample_miss = (dgt == -1).sum(axis=0) / nsites
plt.hist(sample_miss)
sel_sample = sample_miss <= 0.1

samples = samples[sel_sample]
dgt = dgt[:, sel_sample].copy()

# 5. check
nsites, nsamples = dgt.shape
sample_miss = (dgt == -1).sum(axis=0) / nsites
site_miss = (dgt == -1).sum(axis=1) / nsamples
af = (dgt == 1).sum(axis=1) / ((dgt == 0).sum(axis=1) + (dgt == 1).sum(axis=1))

sample_miss.max()
site_miss.max()
af.max()  # =0.992,
af.min()  # =0.007

# 6. final af filter
sel_sites = (af >= 0.01) & (1 - af >= 0.01)
sel_sites.sum()

pos = pos[sel_sites].copy()
chrom = chrom[sel_sites].copy()
alleles = alleles[sel_sites].copy()
dgt = dgt[sel_sites, :].copy()

pos.shape
chrom.shape
dgt.shape
samples.shape

# write to vcf for imputation
import pandas as pd

df = pd.DataFrame(dgt, columns=samples).astype(str)

df = df.replace("-1", "./.")
df = df.replace("1", "1|1")
df = df.replace("0", "0|0")

# TODO: add real REF/ALT
info = pd.DataFrame(
    {
        "#CHROM": chrom,
        "POS": pos,
        "ID": ".",
        "REF": alleles[:, 0],
        "ALT": alleles[:, 1],
        "QUAL": ".",
        "FILTER": "PASS",
        "INFO": ".",
        "FORMAT": "GT",
    }
)
df2 = pd.concat([info, df], axis=1)
header = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=Pf3D7_01_v3,length=640851>
##contig=<ID=Pf3D7_02_v3,length=947102>
##contig=<ID=Pf3D7_03_v3,length=1067971>
##contig=<ID=Pf3D7_04_v3,length=1200490>
##contig=<ID=Pf3D7_05_v3,length=1343557>
##contig=<ID=Pf3D7_06_v3,length=1418242>
##contig=<ID=Pf3D7_07_v3,length=1445207>
##contig=<ID=Pf3D7_08_v3,length=1472805>
##contig=<ID=Pf3D7_09_v3,length=1541735>
##contig=<ID=Pf3D7_10_v3,length=1687656>
##contig=<ID=Pf3D7_11_v3,length=2038340>
##contig=<ID=Pf3D7_12_v3,length=2271494>
##contig=<ID=Pf3D7_13_v3,length=2925236>
##contig=<ID=Pf3D7_14_v3,length=3291936>
##contig=<ID=Pf3D7_API_v3,length=34250>
##contig=<ID=Pf_M76611,length=5967>
"""
df2


with open("pf7_unimp_dom_0_9.vcf", "w") as f:
    f.write(header)
    df2.to_csv(f, sep="\t", index=None, header=True)


cmd = f""" bgzip "pf7_unimp_dom_0_9.vcf" """
run(cmd, shell=True, check=True)
