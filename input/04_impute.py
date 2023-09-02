#! /usr/bin/env python3
from pathlib import Path
from urllib.request import urlretrieve
import re
from subprocess import run

# download beagle
beagle_url = "https://faculty.washington.edu/browning/beagle/beagle.18May20.d20.jar"
beagle_jar = Path("impute/beagle51.jar")
if not beagle_jar.exists():
    beagle_jar.parent.mkdir(parents=True, exist_ok=True)
    urlretrieve(beagle_url, beagle_jar)

# make map files
header_lines = """##contig=<ID=Pf3D7_01_v3,length=640851>
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
##contig=<ID=Pf3D7_14_v3,length=3291936>""".splitlines()


lst = []  # List(chrno, chrname, length)
for line in header_lines:
    pattern = r"ID=(.*?),length=(\d+)"
    match = re.search(pattern, line)
    chrname = match.group(1)
    length = int(match.group(2))
    pattern2 = r"Pf3D7_(\d+)_v3"
    match2 = re.search(pattern2, chrname)
    chrno = int(match2.group(1))
    lst.append((chrno, chrname, length))

bp_per_cm = 15000

with open("impute/map.txt", "w") as f:
    # eg. chrno . cm bp
    for chrno, chrname, chrlen in lst:
        f.write(f"{chrno} . 0.0 1\n")
        f.write(f"{chrno} . {chrlen/bp_per_cm} {chrlen}\n")

# rename chrname names
with open("impute/chrname_map.txt", "w") as f:
    for chrno, chrname, chrlen in lst:
        f.write(f"{chrname}\t{chrno}\n")

# run beagle
input_vcf = "./pf7_unimp_dom.vcf.gz"
output_vcf_prefix = "impute/pf7_dom_imp"
map = "impute/map.txt"

cmd = f"""
    bcftools annotate --rename-chrs impute/chrname_map.txt -Oz -o impute/tmp.vcf.gz {input_vcf}
    java -Xmx10g -jar {beagle_jar} \
            gt=impute/tmp.vcf.gz \
            map={map} \
            out={output_vcf_prefix}_all \
            nthreads=10
    rm impute/tmp.vcf.gz
    bcftools index {output_vcf_prefix}_all.vcf.gz
"""
run(cmd, shell=True, check=True)

# split by chr
for chrno in range(1, 15):
    chr_vcf = output_vcf_prefix + f"_chr{chrno}.vcf.gz"
    cmd = (
        f"""bcftools view -r {chrno} -Oz -o {chr_vcf} {output_vcf_prefix}_all.vcf.gz """
    )
    run(cmd, shell=True, check=True)

for chrno, chrname, chrlen in lst:
    with open(f"impute/map_chr{chrno}.txt", "w") as f:
        # eg. chrno . cm bp
        f.write(f"{chrno} . 0.0 1\n")
        f.write(f"{chrno} . {chrlen/bp_per_cm} {chrlen}\n")
