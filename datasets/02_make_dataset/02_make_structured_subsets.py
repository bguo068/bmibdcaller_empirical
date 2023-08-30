import pandas as pd

unrelated = pd.read_csv("unrelated_samples_10k.txt", names=["Sample"]).Sample
meta = pd.read_csv("../../input/pf7_meta.tsv", sep="\t")

meta = meta[meta.Sample.isin(unrelated)].copy()

meta.Population.value_counts()
# AF-W       5064
# AF-E       1400
# AS-SE-E    1035
# AS-S-FE     994
# AS-SE-W     866
# AF-C        439
# OC-NG       190
# AS-S-E      127
# AF-NE       116
# SA           26

meta[(meta.Year >= 2012) & (meta.Year <= 2013)].groupby(["Population", "Year"])[
    "Sample"
].count().unstack(level=-1).fillna(0).sum(axis=1)[lambda s: s > 30].sum()
# Population
# AF-C       243.0
# AF-E       319.0
# AF-W       994.0
# AS-S-FE     42.0
# AS-SE-E    130.0
# AS-SE-W    180.0
# OC-NG       68.0

sel_year = (meta.Year >= 2012) & (meta.Year <= 2013)
sel_pop = meta.Population.isin("AF-C AF-E AF-W AS-S-FE AS-SE-E AS-SE-W OC-NG".split())
meta2 = meta[sel_pop & sel_year].copy()

# subsampling for population with more than 300 samples
selected_samples = pd.concat(
    [
        meta2[meta2.Population == "AF-C"].Sample,
        meta2[meta2.Population == "AF-E"].Sample.sample(n=300, replace=False),
        meta2[meta2.Population == "AF-W"].Sample.sample(n=300, replace=False),
        meta2[meta2.Population == "AS-S-FE"].Sample,
        meta2[meta2.Population == "AS-SE-E"].Sample,
        meta2[meta2.Population == "AS-SE-W"].Sample,
        meta2[meta2.Population == "OC-NG"].Sample,
    ],
    axis=0,
)

meta3 = meta[meta.Sample.isin(selected_samples)]
meta3.to_csv("meta_structured.txt", index=None, sep="\t")
