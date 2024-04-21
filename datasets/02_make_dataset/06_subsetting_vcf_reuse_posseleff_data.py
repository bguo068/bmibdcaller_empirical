import pandas as pd
from pathlib import Path
from subprocess import run, check_output

meta_all = pd.read_parquet("./posseleff_inputs/meta.parquet")

# ----------------------------------------------
# ------------  Subsetting Vcf files  ----------
# ----------------------------------------------
# ----------------------------------------------

sets = "posseleff_ESEA posseleff_WAF".split()

# iterate over sets
for set0 in sets:
    region = set0.split("_")[1]
    target_dir = Path(f"vcf/{set0}")
    target_dir.mkdir(exist_ok=True, parents=True)

    samples = (
        # Path(f"./posseleff_inputs/unrelated_{region}.txt")
        Path(f"./posseleff_inputs/{region}_samples.txt")
        .read_text()
        .strip()
        .split("\n")
    )

    # write file with a formated file name
    sample_list_fn = target_dir / f"sample_list_{set0}.txt"
    pd.Series(samples).to_csv(sample_list_fn, index=None, header=None)

    # sample map
    sample_map_fn = target_dir / f"sample_name_map_{set0}.txt"
    sample_map = pd.DataFrame(
        {
            "Sample": samples,
            "Id": [f"tsk_{i}" for i, _ in enumerate(samples)],
        }
    )
    sample_map.to_csv(sample_map_fn, sep="\t", index=None, header=None)

    # meta
    meta_fn = f"meta_{set0}.txt"
    meta_all[meta_all.Sample.isin(samples)].to_csv(meta_fn)

    # iterate over chromosomes
    for chrno in range(1, 15):
        vcf_fn = f"./posseleff_inputs/{region}_imputed_fixsamplename.vcf.gz"
        out_vcf_fn = target_dir / f"chr{chrno}.vcf.gz"
        ## 1. subset vcf by samples
        ## 2. remove rare variants or non-segregating variants
        cmd = f"""
        bcftools annotate --rename-chrs posseleff_inputs/chr_name_map.txt {vcf_fn} -Ou \\
          | bcftools view -t {chrno} -Ou \\
          | bcftools view -S {sample_list_fn} -Ou \\
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
