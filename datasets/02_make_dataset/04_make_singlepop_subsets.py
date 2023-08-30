import pandas as pd

meta = pd.read_csv("../../input/pf7_meta.tsv", sep="\t")
unrelated1 = pd.read_csv(
    "unrelated_samples_AS-SE-E.txt", sep="\t", names=["Sample"]
).Sample
unrelated2 = pd.read_csv(
    "unrelated_samples_AF-W.txt", sep="\t", names=["Sample"]
).Sample

# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------

# select samples from AS-SE-E
# Sample by country
meta[meta.Sample.isin(unrelated1)].Country.value_counts().sort_index()
# Cambodia    459
# Laos        286
# Thailand      9
# Vietnam     281

# Sample by year
meta[meta.Sample.isin(unrelated1)].Year.value_counts().sort_index()
# 1993.0      5
# 2005.0      4
# 2007.0     19
# 2008.0     19
# 2009.0     34
# 2010.0    166
# 2011.0    238
# 2012.0    101
# 2013.0     29
# 2014.0     37
# 2015.0     30
# 2016.0     57
# 2017.0    161
# 2018.0    135

# Sample by year and country
sel = (meta.Sample.isin(unrelated1)) & (meta.Year >= 2010)
groups = meta[sel].groupby(["Year", "Country"])
groups["Sample"].count().unstack().fillna(0).astype(int)
# Country  Cambodia  Laos  Thailand  Vietnam
# Year
# 2010.0         90    32         0       44
# 2011.0        119    47         1       71
# 2012.0         52    22         2       25
# 2013.0         26     0         2        1
# 2014.0         23     0         0       14
# 2015.0         10     2         3       15
# 2016.0         44     0         1       12
# 2017.0         24    87         0       50
# 2018.0          4    96         0       35

# subset samples within 2010-2012
is_unrelated = meta.Sample.isin(unrelated1)
is_2010_2012 = (meta.Year >= 2010) & (meta.Year <= 2012)
meta[is_2010_2012 & is_unrelated].shape
# (505, 17)
meta[is_2010_2012 & is_unrelated].to_csv(
    "./meta_singlepop_AS-SE-E_10_12.txt", sep="\t", index=None
)

# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------

# select samples from AF-W
# Sample by country
meta[meta.Sample.isin(unrelated2)].Country.value_counts().sort_index()
# Benin             133
# Burkina Faso       42
# Cameroon          253
# Côte d'Ivoire      67
# Gabon              55
# Gambia            569
# Ghana            2567
# Guinea            135
# Mali              970
# Mauritania         77
# Nigeria            75
# Senegal           121

# Sample by year
meta[meta.Sample.isin(unrelated2)].Year.value_counts().sort_index()
# 1984.0    127
# 1990.0      9
# 1991.0      1
# 2001.0     30
# 2007.0     72
# 2008.0     97
# 2009.0    109
# 2010.0    158
# 2011.0    235
# 2012.0     91
# 2013.0    902
# 2014.0    826
# 2015.0    526
# 2016.0    722
# 2017.0    591
# 2018.0    568

# Sample by year and country
sel = (meta.Sample.isin(unrelated2)) & (meta.Year >= 2009)
groups = meta[sel].groupby(["Year", "Country"])
groups["Sample"].count().unstack().fillna(0).astype(int)

# Country  Benin  Cameroon  Côte d'Ivoire  Gabon  Gambia  Ghana  Guinea  Mali  \
# Year
# 2009.0       0         0              0      0       0    109       0     0
# 2010.0       0         0              0      0       0    158       0     0
# 2011.0       0         0              0      0       0     99     135     0
# 2012.0       0         0              0      0       0     84       0     3
# 2013.0       0       227             49      0      29    239       0   311
# 2014.0      33         0             18     55     193    285       0   121
# 2015.0       0         0              0      0      68    305       0   102
# 2016.0     100         8              0      0      33    242       0   339
# 2017.0       0        18              0      0      24    478       0    22
# 2018.0       0         0              0      0       0    568       0     0

# Country  Mauritania  Nigeria  Senegal
# Year
# 2009.0            0        0        0
# 2010.0            0        0        0
# 2011.0            0        1        0
# 2012.0            0        4        0
# 2013.0            0        0       47
# 2014.0           77       21       23
# 2015.0            0        0       51
# 2016.0            0        0        0
# 2017.0            0       49        0
# 2018.0            0        0        0

# subset samples within 2010-2012
is_unrelated = meta.Sample.isin(unrelated2)
is_2016_2018 = (meta.Year >= 2016) & (meta.Year <= 2018)
is_in_ghana = meta.Country == "Ghana"
meta[is_2016_2018 & is_unrelated & is_in_ghana].shape
# (1288, 17)
meta[is_2016_2018 & is_unrelated & is_in_ghana].to_csv(
    "./meta_singlepop_AF-W_Ghana_16_18.txt", sep="\t", index=None
)
