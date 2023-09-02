import pandas as pd

unrelated = pd.read_csv("unrelated_samples_6k.txt", names=["Sample"]).Sample
meta = pd.read_csv("../../input/pf7_meta.tsv", sep="\t")

meta = meta[meta.Sample.isin(unrelated)].copy()

meta.Population.value_counts()
# AF-W       2523
# AF-E        730
# AS-SE-E     706
# AS-SE-W     639
# AS-S-FE     633
# AF-C        197
# OC-NG       154
# AS-S-E       66
# AF-NE        62
# SA           25

meta[(meta.Year >= 2012) & (meta.Year <= 2013)].groupby(["Population", "Year"])[
    "Sample"
].count().unstack(level=-1).fillna(0).sum(axis=1)[
    lambda s: s > 30
]  # .sum()
# AF-C       103.0
# AF-E       197.0
# AF-W       508.0
# AS-SE-E    100.0
# AS-SE-W    144.0
# OC-NG       57.0

sel_year = (meta.Year >= 2012) & (meta.Year <= 2013)
sel_pop = meta.Population.isin("AF-C AF-E AF-W AS-SE-E AS-SE-W OC-NG".split())
meta2 = meta[sel_pop & sel_year].copy()

# subsampling for population with more than 300 samples
selected_samples = pd.concat(
    [
        meta2[meta2.Population == "AF-C"].Sample,
        meta2[meta2.Population == "AF-E"].Sample,
        meta2[meta2.Population == "AF-W"].Sample.sample(n=300, replace=False),
        # meta2[meta2.Population == "AS-S-FE"].Sample,
        meta2[meta2.Population == "AS-SE-E"].Sample,
        meta2[meta2.Population == "AS-SE-W"].Sample,
        meta2[meta2.Population == "OC-NG"].Sample,
    ],
    axis=0,
)

meta3 = meta[meta.Sample.isin(selected_samples)]
meta3.to_csv("meta_structured.txt", index=None, sep="\t")
