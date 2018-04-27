import mseutils
import pandas as pd
from IPython import embed
from scipy.sparse import csgraph
from collections import Counter

rawdf = pd.read_csv(r"G:\My Drive\_HIFAN_\Data\20140402_RLUS-2048C_iDTs_cppis.csv")
adjmat = mseutils.utils.get_adjacency_matrix_lil(rawdf)
num,labels = csgraph.connected_components(adjmat,directed=False)
rawdf['concomp'] = labels
c = Counter(labels)
biggies = [(l,v) for l,v in c.items() if v > 1]
biggies.sort(key=lambda x:x[1],reverse=True)

biggrp = rawdf[rawdf["concomp"]==biggies[-1][0]]

embed()
print(adjmat.__repr__())