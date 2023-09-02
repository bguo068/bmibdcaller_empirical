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
# Cambodia    324
# Laos        200
# Thailand      3
# Vietnam     179
# Sample by year
meta[meta.Sample.isin(unrelated1)].Year.value_counts().sort_index()
# 1993.0      5
# 2005.0      3
# 2007.0     17
# 2008.0     15
# 2009.0     23
# 2010.0    116
# 2011.0    153
# 2012.0     75
# 2013.0     25
# 2014.0     24
# 2015.0     18
# 2016.0     33
# 2017.0    101
# 2018.0     98

# Sample by year and country
sel = (meta.Sample.isin(unrelated1)) & (meta.Year >= 2010)
groups = meta[sel].groupby(["Year", "Country"])
print(groups["Sample"].count().unstack().fillna(0).astype(int))
# Country  Cambodia  Laos  Thailand  Vietnam
# Year
# 2010.0         61    27         0       28
# 2011.0         82    27         0       44
# 2012.0         43    15         1       16
# 2013.0         23     0         1        1
# 2014.0         16     0         0        8
# 2015.0          3     2         1       12
# 2016.0         24     0         0        9
# 2017.0         16    56         0       29
# 2018.0          2    73         0       23

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
# Benin              85
# Burkina Faso       13
# Cameroon          121
# Côte d'Ivoire      41
# Gabon              33
# Gambia            304
# Ghana            1153
# Guinea             71
# Mali              521
# Mauritania         44
# Nigeria            47
# Senegal            90

# Sample by year
meta[meta.Sample.isin(unrelated2)].Year.value_counts().sort_index()
# 1984.0     28
# 1990.0      2
# 1991.0      1
# 2001.0     15
# 2007.0     36
# 2008.0     50
# 2009.0     41
# 2010.0     80
# 2011.0    111
# 2012.0     43
# 2013.0    465
# 2014.0    474
# 2015.0    281
# 2016.0    354
# 2017.0    267
# 2018.0    275

# Sample by year and country
sel = (meta.Sample.isin(unrelated2)) & (meta.Year >= 2009)
groups = meta[sel].groupby(["Year", "Country"])
groups["Sample"].count().unstack().fillna(0).astype(int)
# print(_)

# Country  Benin  Cameroon  Côte d'Ivoire  Gabon  Gambia  Ghana  Guinea  Mali  \
# Year
# 2009.0       0         0              0      0       0     41       0     0
# 2010.0       0         0              0      0       0     80       0     0
# 2011.0       0         0              0      0       0     39      71     0
# 2012.0       0         0              0      0       0     39       0     2
# 2013.0       0       111             30      0      19     88       0   178
# 2014.0      23         0             11     33     125    135       0    65
# 2015.0       0         0              0      0      46    146       0    57
# 2016.0      62         6              0      0      15     98       0   173
# 2017.0       0         4              0      0      16    212       0    10
# 2018.0       0         0              0      0       0    275       0     0

# Country  Mauritania  Nigeria  Senegal
# Year
# 2009.0            0        0        0
# 2010.0            0        0        0
# 2011.0            0        1        0
# 2012.0            0        2        0
# 2013.0            0        0       39
# 2014.0           44       19       19
# 2015.0            0        0       32
# 2016.0            0        0        0
# 2017.0            0       25        0
# 2018.0            0        0        0


# subset samples within 2016-2018
is_unrelated = meta.Sample.isin(unrelated2)
is_2016_2018 = (meta.Year >= 2016) & (meta.Year <= 2018)
is_in_ghana = meta.Country == "Ghana"
meta[is_2016_2018 & is_unrelated & is_in_ghana].shape
# (1288, 17)
meta[is_2016_2018 & is_unrelated & is_in_ghana].to_csv(
    "./meta_singlepop_AF-W_Ghana_16_18.txt", sep="\t", index=None
)
