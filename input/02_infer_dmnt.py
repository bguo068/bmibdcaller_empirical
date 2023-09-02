import numpy as np

pos = np.load("pos.npy")
chrom = np.load("chr.npy", allow_pickle=True)
samples = np.load("samples.npy", allow_pickle=True)
ad = np.load("/data/bing/pf7_ad.npy")

# ----------------------------------------------
# infer dominant allele
# ---------------------------------------------

# allocate memory
nsites = ad.shape[0]
nsam = ad.shape[1]
gt = np.zeros((nsites, nsam), np.int8)

# init to -1 (missing)
gt[:] = -1

# Calculate dominant genotype based on AD

## Note: due to memory limit, I am doing this by chunks.
chunks = []
i = 0
while i < nsites:
    start = i
    end = i + 20000
    if end > nsites:
        end = nsites
    chunks.append((start, end))
    i = end


for ichunk, (start, end) in enumerate(chunks):
    nsites0 = end - start
    gt0 = np.repeat(-1, nsites0 * nsam).reshape((nsites0, nsam))
    ad0 = ad[start:end, :, 0]
    ad1 = ad[start:end, :, 1]
    total = ad0 + ad1
    ratio = ad0 / total

    sel = (total >= 5) & (ratio >= 0.9)
    gt0[sel] = 0

    sel = (total >= 5) & (ratio <= 0.9)
    gt0[sel] = 1

    # save the chunk to the main gt matrix
    gt[start:end, :] = gt0

    print(f"done {ichunk+1}/{len(chunks)}")

np.save("/data/bing/pf7_gt_dominant.npy", gt)
