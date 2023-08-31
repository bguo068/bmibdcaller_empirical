import pandas as pd
from pathlib import Path
from subprocess import run, check_output


# ----------------------------------------------
# ------------  Subsetting Vcf files  ----------
# ----------------------------------------------
# ----------------------------------------------
def get_in_vcf(chrno: int) -> Path:
    p = Path(f"../../input/impute/pf7_dom_0_9_imp_chr{chrno}.vcf.gz")
    assert p.exists()
    # change suffix from vcf.gz to vcf.gz.csi
    csi = p.with_suffix(".gz.csi")
    # print(csi)
    if not csi.exists():
        run(f"bcftools index {p}", shell=True, check=True)
    assert csi.exists()
    return p


sets = "singlepop_AF-W_Ghana_16_18  singlepop_AS-SE-E_10_12  structured".split()

# iterate over sets
for set0 in sets:
    target_dir = Path(f"vcf/{set0}")
    target_dir.mkdir(exist_ok=True, parents=True)

    sample_list_fn = target_dir / f"sample_list_{set0}.txt"
    meta0 = pd.read_csv(f"meta_{set0}.txt", sep="\t")
    meta0.Sample.to_csv(sample_list_fn, index=None, header=None)
    sample_map_fn = target_dir / f"sample_name_map_{set0}.txt"
    sample_map = pd.DataFrame(
        {"Sample": meta0.Sample, "Id": [f"tsk_{i}" for i in range(meta0.shape[0])]}
    )
    sample_map.to_csv(sample_map_fn, sep="\t", index=None, header=None)

    # iterate over chromosomes
    for chrno in range(1, 15):
        vcf_fn = get_in_vcf(chrno=chrno)
        out_vcf_fn = target_dir / f"chr{chrno}.vcf.gz"
        ## 1. subset vcf by samples
        ## 2. remove rare variants or non-segregating variants
        cmd = f"""
        bcftools view  -S {sample_list_fn} {vcf_fn} -Ou \\
          | bcftools view -q 0.01:minor -Oz \\
          | bcftools reheader -s {sample_map_fn} -o {out_vcf_fn}
        """
        print(cmd)
        run(cmd, shell=True, check=True)

# count number of sites for resulting vcf files (genome-wide)
for set0 in sets:
    target_dir = Path(f"vcf/{set0}")
    target_dir.mkdir(exist_ok=True, parents=True)

    nsites = 0
    for chrno in range(1, 15):
        out_vcf_fn = target_dir / f"chr{chrno}.vcf.gz"
        cmd = "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' " f"{out_vcf_fn} | wc -l"
        nsites += int(check_output(cmd, shell=True, text=True).strip())

    print(f"{set0}: {nsites}")

# ----------------------------------------------
# ------------  counting sites -----------------
# ----------------------------------------------
# ----------------------------------------------
# Number of non-rare sites
# singlepop_AF-W_Ghana_16_18:  15194
# singlepop_AS-SE-E_10_12:     11805
# structured:                  20428
