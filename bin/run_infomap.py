#! /usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np
import pandas as pd
from ibdutils.utils.ibdutils import IBD


def parse_args():
    p = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("--ibd_obj", type=str, required=True)
    # p.add_argument("--npop", type=int, required=True)
    # p.add_argument("--nsam", type=int, required=True)
    p.add_argument("--name_map", type=str, required=True)
    p.add_argument("--genome_set_id", type=int, required=True)
    p.add_argument("--cut_mode", type=str, required=True)
    p.add_argument("--ntrials", type=int, default=1000)
    p.add_argument(
        "--transform", type=str, choices=["square", "cube", "none"], default="square"
    )
    p.add_argument("--ifm_mincm", type=float, default=2)
    p.add_argument("--ifm_mingwcm", type=float, default=5)
    p.add_argument("--ifm_rmchr", type=int, default=0)

    args = p.parse_args()
    if args.transform == "none":
        args.transform = None

    return p.parse_args()


def run(args) -> pd.DataFrame:
    ibd = IBD.pickle_load(args.ibd_obj)

    # tsk_id/PCD113
    name_map = pd.read_csv(args.name_map, sep="\t", names=["RealSample", "Sample"])
    meta = name_map[["Sample", "RealSample"]]
    meta["Sample"] = meta.Sample.str.replace("tsk_", "").astype(int)

    # if ifm_rmchr is a valid chromosome number, remove all IBD from this chromosome
    ibd._df = ibd._df[lambda df: df.Chromosome != args.ifm_rmchr]

    mat = ibd.make_ibd_matrix(min_seg_cm=args.ifm_mincm, min_gw_ibd_cm=args.ifm_mingwcm)
    member_df = ibd.call_infomap_get_member_df(
        mat, meta, trials=args.ntrials, transform=args.transform
    )

    return member_df


if __name__ == "__main__":
    args = parse_args()
    member_df = run(args)

    ofs = f"{args.genome_set_id}_{args.cut_mode}_member.pq"
    member_df.to_parquet(ofs)
    print(member_df)
