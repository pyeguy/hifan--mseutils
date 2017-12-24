from IPython import embed

import mseutils

data_file =  r"G:\My Drive\_HIFAN_\Data\20140402_RLUS-2048C_iDTs_cppis_WPD25_5.csv"

mses = mseutils.file_parsers.load_frag_csv(data_file)

# comb = mseutils.combine_specs(mses)

embed()