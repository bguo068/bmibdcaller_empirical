from ibdutils.utils.ibdutils import IBD
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pandas as pd
import pickle

# resdir of a run of the main pipeline
res_dir = "res/"

simulations = [
    # add a pseudo expected number of peaks
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
    dict(
        genome_set_id=1,
        expected_num_peaks=-1,
        model="mp",
        label="singlepop_AS-SE-E_10_12",
    ),
]

ibdcallers = "hapibd  hmmibd  isorelate  refinedibd  tpbwt".split()

# peaks identified  signal/noise


# ----------------- Signal to noise ratio (per chromosome)   ----------------------
# ----------------- And number of  peaks (expected/observed) ----------------------

# signal/noise ratio results
sn_res = {}
# num of peaks results
np_res = {}
# path to ibd obj
ibdobj_path_res = {}


for simu_ind, simulation in enumerate(simulations):
    # model_ind = 13
    # simulation = simulations[model_ind]
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    model = simulation["model"]

    for caller_ind, ibdcaller in enumerate(ibdcallers):
        # find the right ibdobj file
        ibdobj_fn = ""
        if model == "sp":
            # e.g. 20003_mp_s03/tskibd/ifm_input/20003_orig.ifm.ibdobj.gz
            ibdobj_fn = next(
                Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ibdne_ibd/").glob(
                    "*_orig.ibdne.ibdobj.gz"
                )
            )
        else:
            # eg: 20003_mp_s03/tskibd/ifm_input/20003_orig.ifm.ibdobj.gz
            ibdobj_fn = next(
                Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ifm_input/").glob(
                    "*_orig.ifm.ibdobj.gz"
                )
            )

        # load the data
        tmp = IBD.pickle_load(ibdobj_fn)

        # calculate signal noise ratio
        sn_ratios = []
        for chrno in range(1, 15):
            cov = tmp._cov_df.loc[
                lambda df: df.Chromosome == chrno, ["Start", "Coverage"]
            ].set_index("Start")
            q5, q95 = np.quantile(cov, [0.05, 0.95])
            noise = cov[(cov >= q5) & (cov <= q95)].mean().to_numpy()[0]
            signal = cov.max().to_numpy()[0]
            sn_ratio = signal / noise
            sn_ratios.append(sn_ratio)
        sn_res[(label, ibdcaller)] = sn_ratios

        num_peaks = tmp._peaks_df.shape[0]
        expected = simulation["expected_num_peaks"]

        np_res[(label, ibdcaller)] = [num_peaks, expected]

        ibdobj_path_res[(label, ibdcaller)] = ibdobj_fn


# ----------------- Ne difference  ------------------------------------------
ne_res = {}

for simu_ind, simulation in enumerate(simulations):
    # simu_ind = 3
    # simulation = simulations[simu_ind]
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    model = simulation["model"]

    # Ne is estimated only for SP model
    if model != "sp":
        continue

    for caller_ind, ibdcaller in enumerate(ibdcallers):
        # caller_ind = 0
        # ibdcaller = ibdcallers[caller_ind]
        # print(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ne_output")
        fn1 = next(
            Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ne_output").glob(
                "*_orig.ne"
            )
        )
        df = pd.read_csv(fn1, sep="\t")
        df.columns = ["GEN", "NE", "L95", "U95"]
        ne_res[(label, ibdcaller, "orig")] = df

        fn2 = next(
            Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ne_output").glob(
                "*_rmpeaks.ne"
            )
        )
        df = pd.read_csv(fn2, sep="\t")
        df.columns = ["GEN", "NE", "L95", "U95"]
        ne_res[(label, ibdcaller, "rmpeaks")] = df


# Infomap difference
ifm_res = {}
for simu_ind, simulation in enumerate(simulations):
    # simu_ind = 3
    # simulation = simulations[simu_ind]
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    model = simulation["model"]

    # Infomap assignment is estimated only for MP model
    if model != "mp":
        continue

    for caller_ind, ibdcaller in enumerate(ibdcallers):
        # caller_ind = 0
        # ibdcaller = ibdcallers[caller_ind]
        # e.g. 20003_mp_s03/tskibd/ifm_output/20003_orig_member.pq
        # print(Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ifm_output"))
        fn1 = next(
            Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ifm_output").glob(
                "*_orig_member.pq"
            )
        )
        df = pd.read_parquet(fn1)
        ifm_res[(label, ibdcaller, "orig")] = df

        fn2 = next(
            Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ifm_output").glob(
                "*_rmpeaks_member.pq"
            )
        )
        df = pd.read_parquet(fn1)
        ifm_res[(label, ibdcaller, "rmpeaks")] = df


# IBD comparison
raw_ibd_file_res = {}

for simu_ind, simulation in enumerate(simulations):
    # simu_ind = 3
    # simulation = simulations[simu_ind]
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    model = label.split("_")[0]

    for caller_ind, ibdcaller in enumerate(ibdcallers):
        # caller_ind = 0
        # ibdcaller = ibdcallers[caller_ind]
        # e.g.  sp_s03/ibd/tskibd/10003_10_tskibd.ibd
        ibd_fn_lst = []
        for chrno in range(1, 15):
            fn = f"{res_dir}/{label}/ibd/{ibdcaller}/{genome_set_id}_{chrno}_{ibdcaller}.ibd"
            # fix inconsistent naming pattern
            if ibdcaller == "tpbwt":
                fn = f"{res_dir}/{label}/ibd/{ibdcaller}ibd/{genome_set_id}_{chrno}_{ibdcaller}ibd.ibd"
            assert Path(fn).exists()
            ibd_fn_lst.append(chrno)
        assert Path(fn).exists()
        raw_ibd_file_res[(label, ibdcaller)] = ibd_fn_lst

# write file
outdir = Path("./analysis_res/info")
outdir.mkdir(parents=True, exist_ok=True)


def pickle_to_file(obj, fn):
    with open(fn, "wb") as f:
        pickle.dump(obj, f)


pickle_to_file(ibdobj_path_res, f"{outdir}/ibdobj_path.pkl")
pickle_to_file(raw_ibd_file_res, f"{outdir}/raw_ibd_path.pkl")
pickle_to_file(sn_res, f"{outdir}/signal_noise_ratio.pkl")
pickle_to_file(np_res, f"{outdir}/num_peaks.pkl")
pickle_to_file(ne_res, f"{outdir}/ne_res.pkl")
pickle_to_file(ifm_res, f"{outdir}/infomap_res.pkl")
