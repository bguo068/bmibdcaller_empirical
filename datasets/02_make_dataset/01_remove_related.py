from ibdutils.utils.ibdutils import IBD, Genome

ibd_fn_lst = [
    f"../01_call_hapibd/res/alldom09imp/ibd/hapibd/0_{chrno}_hapibd.ibd"
    for chrno in range(1, 15)
]

genome = Genome.get_genome("Pf3D7")
ibd = IBD(genome, label="alldom09imp")
ibd.read_ibd(ibd_fn_lst)

M = ibd.make_ibd_matrix()
unrelated = ibd.get_unrelated_samples(M)

unrelated.to_csv("unrelated_samples_10k.txt", index=None, header=None)
