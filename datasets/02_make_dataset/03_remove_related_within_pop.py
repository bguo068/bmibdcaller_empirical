from ibdutils.utils.ibdutils import IBD, Genome
import pandas as pd

# load meta for all samples
meta = pd.read_csv("../../input/pf7_meta.tsv", sep="\t")

# make a list of ibd files
ibd_fn_lst = [
    f"../01_call_hapibd/res/alldom09imp/ibd/hapibd/0_{chrno}_hapibd.ibd"
    for chrno in range(1, 15)
]

# construct IBD object
genome = Genome.get_genome("Pf3D7")
ibd = IBD(genome, label="alldom09imp")
ibd.read_ibd(ibd_fn_lst)

# show the number of samples in each population
meta.Population.value_counts()

# subset samples by a specific population


def subset_for_as_se_e(meta, ibd):
    meta1 = meta[meta["Population"] == "AS-SE-E"]
    pop_samples = meta1.Sample
    ibd1 = ibd.duplicate(new_label="alldom09imp_AS-SE-E")
    ibd1.subset_ibd_by_samples(pop_samples)

    M1 = ibd1.make_ibd_matrix()
    unrelated1 = ibd1.get_unrelated_samples(M1)

    unrelated1.to_csv("unrelated_samples_AS-SE-E.txt", index=None, header=None)


subset_for_as_se_e(meta, ibd)


def subset_for_af_w():
    meta2 = meta[meta["Population"] == "AF-W"]
    pop_samples = meta2.Sample
    ibd2 = ibd.duplicate(new_label="alldom09imp_AF-W")
    ibd2.subset_ibd_by_samples(pop_samples)

    M2 = ibd2.make_ibd_matrix()
    unrelated2 = ibd2.get_unrelated_samples(M2)

    unrelated2.to_csv("unrelated_samples_AF-W.txt", index=None, header=None)

    meta[meta.Sample.isin(unrelated2)]

    meta[meta.Sample.isin(unrelated2)].Country.value_counts()
    meta[meta.Sample.isin(unrelated2)].Year.value_counts().sort_index()
    meta[(meta.Sample.isin(unrelated2)) & (meta.Year > 2009)].groupby(
        ["Country", "Year"]
    )["Sample"].count().unstack().fillna(0).astype(int)

    meta[(meta.Sample.isin(unrelated2)) & (meta.Year >= 2016) & (meta.Year <= 2018)]


subset_for_af_w(meta, ibd)
