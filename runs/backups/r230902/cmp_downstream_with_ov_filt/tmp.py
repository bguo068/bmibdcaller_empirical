from ibdutils.utils.ibdutils import IBD
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pathlib import Path
import numpy as np
import pandas as pd
import pickle
import tomli_w
from subprocess import run
from shutil import rmtree
import igraph as ig

indir = Path("./analysis_res/info")
indir.mkdir(parents=True, exist_ok=True)


def pickle_load_file(fn):
    with open(fn, "rb") as f:
        return pickle.load(f)


# get Sample, population from the general meta file
meta = pd.read_csv("../../../input/pf7_meta.tsv", sep="\t")
ibdobj_path_res = pickle_load_file(f"{indir}/ibdobj_path.pkl")
sn_res = pickle_load_file(f"{indir}/signal_noise_ratio.pkl")
np_res = pickle_load_file(f"{indir}/num_peaks.pkl")
ne_res = pickle_load_file(f"{indir}/ne_res.pkl")
ifm_res = pickle_load_file(f"{indir}/infomap_res.pkl")

simulations = [
    # add a pseudo expected number of peaks the following 3 lines
    dict(
        genome_set_id=0,
        expected_num_peaks=-1,
        model="sp",
        label="singlepop_AF-W_Ghana_16_18",
    ),
    dict(
        genome_set_id=1,
        expected_num_peaks=-1,
        model="sp",
        label="singlepop_AS-SE-E_10_12",
    ),
    dict(genome_set_id=100, expected_num_peaks=-1, model="mp", label="structured"),
]

ibdcallers = "hapibd  hmmibd  isorelate  refinedibd  tpbwt".split()
ibdcallers = "refinedibd hmmibd".split()


def get_sample_count(label):
    p = Path(
        f"../../../datasets/02_make_dataset/vcf/{label}/sample_name_map_{label}.txt"
    )
    assert p.exists()
    df = pd.read_csv(p, sep="\t", header=None)
    return df.shape[0]


# outdir = Path("./analysis_res")
# outdir.mkdir(exist_ok=True, parents=True)

# ------------------------- plot coverage and peaks --------------------------

for simulation in simulations:
    label = simulation["label"]
    fig, axes = plt.subplots(
        nrows=5,
        ncols=1,
        figsize=(10, 8),
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )
    for icaller, caller in enumerate(ibdcallers):
        ax = axes[icaller]
        if icaller == 0:
            ax.set_title(label)
        ibd = IBD.pickle_load(ibdobj_path_res[label, caller])
        ibd.plot_coverage(ax, label=caller, which="xirsfilt", plot_proportions=False)
    # fig.savefig(f"{outdir}/coverage_{label}.png", dpi=600)

ibd = IBD.pickle_load(ibdobj_path_res["singlepop_AF-W_Ghana_16_18", "refinedibd"])
cm = ibd.calc_ibd_length_in_cm()

# -------------------------------------------------------------

df = pd.DataFrame(
    [
        {
            "Label": label,
            "Caller": ibdcaller,
            "SNRatio": np.mean(r_lst),
            "Std": np.std(r_lst),
        }
        for (label, ibdcaller), r_lst in sn_res.items()
    ]
)

fig = plt.figure(constrained_layout=True, figsize=(8, 8))
for ilabel, simulation in enumerate(simulations):
    label = simulation["label"]
    df_sel = df[df.Label == label]
    ax = fig.add_subplot(3, 2, ilabel + 1)
    ax.bar(df_sel.Caller, df_sel.SNRatio, yerr=df_sel.Std, capsize=5)
    ax.set_ylabel("signal to noise ratio")
    ax.set_title(f"simulation: {label}")
    ax.set_xticklabels(
        ax.get_xticklabels(), rotation=45, rotation_mode="anchor", va="top", ha="right"
    )

fig.savefig(f"{outdir}/signal_noise.png", dpi=600)

# ------------------------------------------------------------------------

df = dict(Label=[], Caller=[], Obs=[], Exp=[])
for (label, caller), (obs, exp) in np_res.items():
    df["Label"].append(label)
    df["Caller"].append(caller)
    df["Obs"].append(obs)
    df["Exp"].append(exp)
df = pd.DataFrame(df)

fig = plt.figure(constrained_layout=True, figsize=(8, 8))
for ilabel, simulation in enumerate(simulations):
    label = simulation["label"]
    # label = "mp_s01"
    df_sel = df[df.Label == label]
    ax = fig.add_subplot(4, 2, ilabel + 1)
    ax.bar(df_sel.Caller, df_sel.Obs)
    ax.set_title(f"simulation: {label}")
    ax.set_xticklabels(
        ax.get_xticklabels(), rotation=45, rotation_mode="anchor", va="top", ha="right"
    )
    ax.set_ylabel("no. of peaks observed")
outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/num_peaks_observed.png", dpi=600)

# --------------------------- NE -------------------------------------


fig = plt.figure(figsize=(15, 5), constrained_layout=True)
# only the first sets are treated as single population
for ilabel, simulation in enumerate(simulations[:2]):
    target_label = simulation["label"]
    # target_label = "sp_neu"
    target_rm_model = "orig"

    param_df = None
    df_res = {}
    # filter ne dataframes for a given
    for (label, caller, rm_mode), df in ne_res.items():
        if label != target_label:
            continue
        if rm_mode != target_rm_model:
            continue
        else:
            df_res[caller] = df
            print(label, caller, rm_mode)
    df_res.keys()

    ncallers = len(ibdcallers)
    for icaller, caller in enumerate(ibdcallers):
        ax = fig.add_subplot(2, ncallers, ilabel * ncallers + icaller + 1)
        df_sel = df_res[caller]
        # ax.plot(param_df.GEN, param_df.NE, label="Truth", color="k", linestyle="--")
        ax.plot(df_sel.GEN, df_sel.NE, label=caller, color="r")
        ax.fill_between(df_sel.GEN, y1=df_sel.L95, y2=df_sel.U95, fc="r", alpha=0.3)
        ax.set_xlim(0, 100)
        ax.set_ylim(1e2, 1e6)
        ax.set_title(caller)
        ax.legend()
        ax.set_yscale("log")
        if icaller == 0:
            ax.set_ylabel(f"{target_label.upper()}\n\nNe")
        if ilabel == 1:
            ax.set_xlabel("generations ago")
fig.suptitle("rm_mode = orig")

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/ne.png", dpi=600)


# --------------------------- infomap -------------------------------------
def transform_ifm_df(df):
    """transform infomap df to a 5x5 matrix

    args:
    -----
    df: infomap df (contains Sample, Rank, Population columns)

    columns: community label
    rows: true population label
    cells: number of samples with a given true population label and a given
    community label

    """
    npop = df.Population.unique().size
    df = (
        df.groupby(["Population", "Rank"])["Sample"]
        .count()
        .unstack()
        .fillna(0)
        .astype(int)
        .iloc[:, :npop]
        .copy()
    )
    # make up columns if necessary
    if df.shape[1] < npop:
        for i in range(df.shape[1], npop):
            df[f"c{i}"] = 0
    df = df.sort_values(by=df.index.to_list(), axis=1, ascending=False)
    df.columns = [f"c{i+1}" for i in range(0, npop)]
    # df.index = [f"p{i+1}" for i in range(0, npop)]
    return df


def get_adj_rank(df):
    """get adjusted rand index from infomap df"""
    pop = pd.Categorical(df.Population)
    pop = pop.codes
    adj_rand = ig.compare_communities(df.Rank, pop, method="adjusted_rand")
    return adj_rand


def merge_with_meta(df, meta):
    """merge df with meta"""
    df = df.copy()
    df["Sample"] = df["RealSample"]
    df = df.merge(meta[["Sample", "Population"]], on="Sample", how="left")
    return df


fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16, 3), constrained_layout=True)
for ilabel, label in enumerate(["structured"]):
    adj_rank_lst = []
    for icaller, caller in enumerate(ibdcallers):
        ifm = ifm_res[(label, caller, "orig")]
        ifm = merge_with_meta(ifm, meta)
        npop = ifm.Population.unique().size
        adj_rank = get_adj_rank(ifm)
        adj_rank_lst.append(adj_rank)
        ax = axes[icaller]
        df = transform_ifm_df(ifm)
        ax.imshow(
            df, cmap="Blues", vmin=0, vmax=200, interpolation="none", aspect="auto"
        )
        ax.set_xticks(np.arange(0, npop))
        ax.set_xticklabels(df.columns)
        ax.set_yticks(np.arange(0, npop))
        if ilabel == 0:
            ax.set_title(caller)
        if ilabel == 3:
            ax.set_xlabel("community label")
        if icaller == 0:
            ax.set_ylabel(f"{label}\n\ntrue pop label")
            ax.set_yticklabels(df.index)
    ax = axes[len(ibdcallers)]
    ax.bar(np.arange(len(ibdcallers)), adj_rank_lst)
    ax.set_xticks(np.arange(len(ibdcallers)))
    ax.set_xticklabels(ibdcallers, rotation=45, rotation_mode="anchor", ha="right")
    ax.set_ylabel("adjusted rand index")
    ax.set_ylim(0, 1)
    if ilabel == 0:
        ax.set_title("adjusted rand index")
fig.savefig(f"{outdir}/infomap.png", dpi=600)


# ------------------------- IBD comparision -----------------------------------
chrlens = [
    640851,
    947102,
    1067971,
    1200490,
    1343557,
    1418242,
    1445207,
    1472805,
    1541735,
    1687656,
    2038340,
    2271494,
    2925236,
    3291936,
]

# prepare map files
bp_per_cm = 15000
r = 0.01 / bp_per_cm
for chrno in range(1, 15):
    len_bp = chrlens[chrno - 1]
    len_cm = len_bp / bp_per_cm
    map_fn = f"tmp/map/{chrno}.map"
    Path(map_fn).parent.mkdir(parents=True, exist_ok=True)
    with open(map_fn, "w") as f:
        f.write(f"{chrno} . 0.0 1\n")
        f.write(f"{chrno} . {len_cm} {len_bp}\n")

# prepare gnome.toml file
genome = {}
# genome file path
# "./r0003_genome.toml"
genome["name"] = f"r{int(r*1e9):04}"
genome["chromsize"] = [
    int(l * 1.1) for l in chrlens
]  # make it larger to avoid issues due rounding errors
genome["chromnames"] = [f"{i}" for i in range(1, 15)]
genome["gmaps"] = [f"tmp/map/{chrno}.map" for chrno in range(1, 15)]
genome["idx"] = {f"{i}": i - 1 for i in range(1, 15)}

genome_fn = f"./tmp/genome.toml"
with open(genome_fn, "wb") as f:
    tomli_w.dump(genome, f)


root_dir = "."

for simulation in simulations:
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]

    # prepare sample file
    nsam = get_sample_count(label)
    with open("tmp/samples.txt", "w") as f:
        for i in range(0, nsam):
            f.write(f"{i}\n")

    # use hmmibd as a proxy for tskibd
    dir = Path(f"tmp/{label}/hmmibd")
    if dir.exists():
        rmtree(dir)
    # link true ibd
    for chrno in range(1, 15):
        # res/mp_s00/ibd/hmmibd/20000_11_hmmibd.ibd
        src_fn = f"{root_dir}/res/{label}/ibd/hmmibd/{genome_set_id}_{chrno}_hmmibd.ibd"
        dst_fn = f"tmp/{label}/hmmibd/{chrno}.ibd"
        Path(dst_fn).parent.mkdir(parents=True, exist_ok=True)
        Path(dst_fn).hardlink_to(src_fn)

    for caller in ibdcallers:
        if caller == "hmmibd":
            continue
        if caller == "tpbwt":
            caller_ = "tpbwtibd"
        else:
            caller_ = caller

        dir = Path(f"tmp/{label}/{caller}/")
        if dir.exists():
            rmtree(dir)
        # link inferred ibd
        for chrno in range(1, 15):
            # res/mp_s00/ibd/hmmibd/20000_11_hmmibd.ibd
            src_fn = f"{root_dir}/res/{label}/ibd/{caller_}/{genome_set_id}_{chrno}_{caller_}.ibd"
            dst_fn = f"tmp/{label}/{caller}/{chrno}.ibd"
            Path(dst_fn).parent.mkdir(parents=True, exist_ok=True)
            Path(dst_fn).hardlink_to(src_fn)

        out_prefix = f"tmp/cmp/{label}/hmmibd_vs_{caller}.tsv"
        Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

        # call ishare/ibdutils
        trueibd_dir = f"tmp/{label}/hmmibd/"
        infribd_dir = f"tmp/{label}/{caller}/"
        cmd = f"""
        /data/bing/ishare/target/release/ibdutils compare \\
            -g tmp/genome.toml \\
            -s tmp/samples.txt \\
            -S tmp/samples.txt \\
            -f tskibd \\
            -F tskibd \\
            -i {trueibd_dir} \\
            -I {infribd_dir} \\
            -o {out_prefix}
        """
        ret = run(cmd, shell=True, text=True)
        if ret.returncode != 0:
            print("error in running ishare/ibdutils")
            print(f"command: \n {cmd}")
            print(f"\n Error:")
            print(ret.stderr)
            raise ("Error")


# read stats
df_lst = []
for simulation in simulations:
    label = simulation["label"]
    for caller in ibdcallers:
        if caller == "hmmibd":
            continue
        csv_fn = f"tmp/cmp/{label}/hmmibd_vs_{caller}.tsv"
        df = pd.read_csv(csv_fn, sep=",")
        df["Label"] = label
        df["Caller"] = caller
        df_lst.append(df)
df = pd.concat(df_lst, axis=0)

nrows = 3
ncols = 4
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(10, 6),
    constrained_layout=True,
    sharex=True,
    sharey=True,
)
for row, simulation in enumerate(simulations):
    label = simulation["label"]
    for col, caller in enumerate([c for c in ibdcallers if c != "hmmibd"]):
        df_sel = df[(df.Label == label) & (df.Caller == caller)].copy()
        df_sel["FN"] = 1 - df_sel.RateAOverlapByB
        df_sel["FP"] = 1 - df_sel.RateBOverlapByA
        df_sel = df_sel.set_index("LenBinStart")
        #
        ax = axes[row, col]
        # ax = fig.add_subplot(nrows, ncols, 1 + row * ncols + col)
        bin_map = {
            "3": "[3-4)",
            "4": "[4-6)",
            "6": "[6-10)",
            "10": "[10-18)",
            "18": "[18-Inf)",
            "genome_wide": "gw",
        }
        for i, b in enumerate(["3", "4", "6", "10", "18", "genome_wide"]):
            marker = list("osd*hx")[i]
            color = "red blue green purple orange brown cyan".split()[i]
            x, y = df_sel.loc[b, "FN"], df_sel.loc[b, "FP"]
            ax.plot(
                [x],
                [y],
                label=bin_map[b],
                marker=marker,
                color=color,
                ms=5,
                linestyle="none",
            )
        if col == ncols - 1:
            ax.legend(
                ncol=2,
                borderpad=0.0,
                labelspacing=0.1,
                handletextpad=0.4,
                handlelength=0.1,
                bbox_to_anchor=(1, 0.5),
            )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        if row == nrows - 1:
            ax.set_xlabel("FN")
        if col == 0:
            ax.set_ylabel(f"{label}\n\nFP")
        if row == 0:
            ax.set_title(caller)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
fig.savefig(f"{outdir}/ibd_cmp_overlap.png", dpi=600)

# ---

totibd_res = {}

for simulation in simulations:
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    #
    for caller in ibdcallers:
        if caller == "hmmibd":
            continue
        out_prefix = f"tmp/cmp/{label}/hmmibd_vs_{caller}.totibd"
        df = pd.read_parquet(out_prefix)

        totibd_res[(label, caller)] = df


nrows = 3
ncols = 4
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(10, 8),
    constrained_layout=True,
    sharex=True,
    sharey=True,
)
for row, simulation in enumerate(simulations):
    label = simulation["label"]
    for col, caller in enumerate([c for c in ibdcallers if c != "hmmibd"]):
        df_sel = totibd_res[(label, caller)]
        df_sel.columns = ["True", "Infer"]
        df_sel = df_sel.loc[lambda df: ~((df["True"] == 0) & (df.Infer == 0)), :]
        #
        ax = axes[row, col]
        ax.plot(df_sel["True"], df_sel["Infer"], ".", alpha=0.3)
        ax.plot([2, 1400], [2, 1400])
        ax.set_yscale("log")
        ax.set_xscale("log")
        # ax.legend(
        #     ncol=2, borderpad=0.0, labelspacing=0.1, handletextpad=0.4, handlelength=0.1
        # )
        if row == nrows - 1:
            ax.set_xlabel("true total IBD")
        if col == 0:
            ax.set_ylabel(f"{label}\n\ninferred total IBD")
        if row == 0:
            ax.set_title(caller)
        # ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        # ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/ibd_cmp_totibd.png", dpi=600)

# -------------------------- plot peaks ----------------
label = "singlepop_AF-W_Ghana_16_18"
caller = "hapibd"
ibd = IBD.pickle_load(ibdobj_path_res[(label, caller)])
fig, ax = plt.subplots(figsize=(12, 3), constrained_layout=True)
ibd.plot_coverage(ax=ax)


# -------------------------- binned population-level total IBD ----------------
def get_bin_pop_total_ibd(label, caller):
    ibd = IBD.pickle_load(ibdobj_path_res[(label, caller)])
    cm = ibd.calc_ibd_length_in_cm()
    bins = np.arange(2, 100, 0.05)
    bin_center = (bins[:-1] + bins[1:]) / 2
    bin_count, _ = np.histogram(cm, bins=bins)
    return (bin_center, bin_center * bin_count)


nrows = 3
ncols = 4
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(10, 8),
    constrained_layout=True,
    sharex=True,
    sharey=True,
)
for row, simulation in enumerate(simulations):
    label = simulation["label"]
    bin_center, x0 = get_bin_pop_total_ibd(label, "hmmibd")
    for col, caller in enumerate([c for c in ibdcallers if c != "hmmibd"]):
        # true
        x = x0.copy()
        # infer
        _, y = get_bin_pop_total_ibd(label, caller)
        sel1 = (x == 0) & (y == 0)
        sel2 = (x == 0) & (y != 0)
        sel3 = (x != 0) & (y == 0)
        sel4 = (x != 0) & (y != 0)
        x[sel2] = 1
        y[sel3] = 1
        ax = axes[row, col]
        ax.plot(bin_center[sel4], y[sel4] / x[sel4], ".", color="blue", alpha=0.3)
        ax.plot(
            bin_center[sel2 | sel3],
            y[sel2 | sel3] / x[sel2 | sel3],
            ".",
            color="red",
            alpha=0.3,
        )
        ax.axhline(y=1.0, color="r", linestyle="--")
        ax.set_yscale("log")
        ax.set_xscale("log")
        if row == nrows - 1:
            ax.set_xlabel("true total IBD")
        if col == 0:
            ax.set_ylabel(f"{label}\n\ninferred total IBD")
        if row == 0:
            ax.set_title(caller)
        # ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        # ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
fig.suptitle("Population-level total length of IBD of 0.05cM length windows")

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/ibd_cmp_pop_level_totibd.png", dpi=600)
