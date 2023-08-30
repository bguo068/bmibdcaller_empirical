import pandas as pd

tskibd = pd.read_csv("./relevant_files/10000_1_tskibd.ibd", sep="\t")
hmmo = pd.read_csv("./relevant_files/hmmibd_out.hmm.txt", sep="\t")
hmmf = pd.read_csv("./relevant_files/hmmibd_out.hmm_fract.txt", sep="\t")

hmmo = hmmo[lambda df: (df["end"] - df["start"] > 15000 * 2) & (df["different"] == 0)]
hmmo = hmmo[["sample1", "sample2", "start", "end"]].copy()

hmmo.columns = ["Id1", "Id2", "Start", "End"]
x = hmmo.Id1.str.replace("tsk_", "").astype(int)
y = hmmo.Id2.str.replace("tsk_", "").astype(int)
xy = pd.concat([x, y], axis=1)
id1 = xy.max(axis=1)
id2 = xy.min(axis=1)
hmmo["Id1"] = id1
hmmo["Id2"] = id2

hmmf = hmmf[["sample1", "sample2", "N_generation"]].copy()
hmmf.columns = ["Id1", "Id2", "Tmrca"]
x = hmmf.Id1.str.replace("tsk_", "").astype(int)
y = hmmf.Id2.str.replace("tsk_", "").astype(int)
xy = pd.concat([x, y], axis=1)
id1 = xy.max(axis=1)
id2 = xy.min(axis=1)
hmmf["Id1"] = id1
hmmf["Id2"] = id2

hmmibd = hmmo.merge(hmmf, on=["Id1", "Id2"], how="left")
hmmibd

hmmibd[lambda df: df.Tmrca < 1.5].sort_values(["Id1", "Id2"])
tskibd[lambda df: df.Tmrca < 1.5].sort_values(["Id1", "Id2"])
a = tskibd[lambda df: df.Tmrca < 1.5]
b = hmmibd[lambda df: df.Tmrca < 3]

amb = a.merge(b, on=["Id1", "Id2"], how="left")

# hmmibd capure 72.5 of the IBD segments with tmrca<1.5 given length >90
(
    amb[lambda df: df.End_x - df.Start_x > 15000 * 90]
    .assign(Hit=lambda df: df.End_y.notna())
    .groupby(["Id1", "Id2"])["Hit"]
    .sum()
    > 0
).mean()
# hmmibd capure 40.1 of the IBD segments with tmrca<1.5 given length >80
(
    amb[lambda df: df.End_x - df.Start_x > 15000 * 80]
    .assign(Hit=lambda df: df.End_y.notna())
    .groupby(["Id1", "Id2"])["Hit"]
    .sum()
    > 0
).mean()

import matplotlib.pyplot as plt

# It seems that there are quite of fraction of the true IBD segmetns with
# tmrca<1.5 are very long.  Given that the error rate of very long IBD segments
# is relatively low.  We many be able use hmmIBD as a alternative to tskit IBD
# for the purpose of detecting IBD segments with tmrca<1.5.
a_cm = (a.End - a.Start) / 15000
plt.hist(a_cm)


# 351 segments with length > 50cM
tskibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 50]
# 202 segments with length > 60cM
tskibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 60]
# 141 segments with length > 70cM
tskibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 70]
# 101 segments with length > 80cM
tskibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 80]
# 54 segments with length > 90cM
tskibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 90].shape

# 349 segments with length > 50cM
hmmibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 50]
# 208 segments with length > 60cM
hmmibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 60]
# 143 segments with length > 70cM
hmmibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 70]
# 103 segments with length > 80cM
hmmibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 80]
# 57 segments with length > 90cM
hmmibd.assign(Cm=lambda df: (df.End - df.Start) / 15000)[lambda df: df.Cm > 90].shape

# two possible soluations:
# 1. use hmmIBD to detect IBD segments with tmrca<1.5
# 2. directly exclude IBD segments that are extremly long like half of the chromosome
