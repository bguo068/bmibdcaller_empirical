It seems that removing inferred IBD that overlapped with true IBD with trmca < 1.5
can largely address the oscillation issue observed in IBDNe results.

However, with empirical dataset, genealogical trees are not available, so we cannot. 
We might be able to get the trmca information from hmmibd results.

This folder aims to check how hmmibd `N_generation` (equivalent to trmca) is reliable
compared with tmrca from tskibd.