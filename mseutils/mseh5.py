'''
Work in progress...

This should be a mapping to HDF5 file format which I think is a natural fit for this MS data..
Might be overkill..

Basic Structure - 
root /
    MSeData / 
        samp1/
            
        samp2
        samp3
        ... 


'''

MS2_ARR_LEN = 200
MGF_FILE_STR_LEN = 5000
SRC_FRG_ARR_LEN = 50 

import tables as tbls

from . import *

from collections import namedtuple

h5_ret_tup = namedtuple('h5_ret_tup','h5 group parent_table srcfrg_table')

def create_h5_file(fname,title=''):
    '''
    Creates a h5 file on disk. Will overwrite old file with same name.

    Args:
        fname (str) = the filename that will be used as the h5 file
        title (str) = optional title to include in the file
    Returns:
        h5_ret_tup : a named tuple which has the hdf5 object (h5) group, 
            parent_table, and source
    '''
    h5 = tbls.open_file(fname,mode='w',title=title)
    group = h5.create_group("/","msedata","MSe Data")
    parent_table = h5.create_table(group, "mse_specs",MSE,"MSe Spectra")
    srcfrg_tbl = h5.create_table(group,'source_frags',MSE,"Source Fragments")
    return h5_ret_tup(h5,group,parent_table,srcfrg_tbl)

def open_h5_file(fname,mode,group_name,parent_table_name,srcfrg_tbl_name):
    '''return a h5_ret_tup'''
    h5 = tbls.open_file(fname,mode)
    group = h5.get_node("/{}".format(group_name))
    parent_table = h5.get_node("/{}/{}".format(group_name,parent_table_name))
    srcfrg_tbl = h5.get_node("/{}/{}".format(group_name,srcfrg_tbl_name))
    return h5_ret_tup(h5,group,parent_table,srcfrg_tbl)

def add_mse(mse,idx,parent_t,srcfrg_t):
    '''
    lots of shitty encoding decoding

    '''
    r = parent_t.row
    r['sampid'] = mse.sampid
    r['idx'] = idx
    r['rt'] = mse.rt.val
    r['ccs'] = mse.ccs.val
    r['mz'] = mse.mz.mz
    r['ppm'] = mse.mz.ppm
    r['z'] = mse.mz.z
    r['n'] = mse.n
    r['intensity'] = mse.i
    ms2arr = [[mz.mz,i] for mz,i in mse.ms2_data.items()]
    ms2arr.extend([[0,0] for _ in range(MS2_ARR_LEN-len(ms2arr))])
    r['ms2_data'] = ms2arr
    r['mgf_files'] = "|".join(mse.mgf_files)
    srcfrgarr = [idx+i+1 for i in range(len(mse.src_frags))]
    srcfrgarr.extend([0 for _ in range(SRC_FRG_ARR_LEN-len(srcfrgarr))])
    r['src_frag_ids'] = srcfrgarr
    r.append()

    # should this be recursive?
    for srcfrg in mse.src_frags:
        r = srcfrg_t.row
        idx +=1
        r['sampid'] = mse.sampid
        r['idx'] = idx
        r['rt'] = srcfrg.rt.val
        r['ccs'] = srcfrg.ccs.val
        r['mz'] = srcfrg.mz.mz
        r['ppm'] = srcfrg.mz.ppm
        r['z'] = srcfrg.mz.z
        r['n'] = srcfrg.n
        r['intensity'] = srcfrg.i
        ms2arr = [[mz.mz,i] for mz,i in mse.ms2_data.items()]
        ms2arr.extend([[0,0] for _ in range(MS2_ARR_LEN-len(ms2arr))])
        r['ms2_data'] = ms2arr
        r['mgf_files'] = "|".join(srcfrg.mgf_files)
        srcfrgarr = [idx+i+1 for i in range(len(mse.src_frags))]
        srcfrgarr.extend([0 for _ in range(SRC_FRG_ARR_LEN-len(srcfrgarr))])
        r['src_frag_ids'] = srcfrgarr        
        r.append()
    return idx

def add_mses(mses,h5t,idx=0,sampid='n/a'):
    parent_t = h5t.parent_table
    srcfrg_t = h5t.srcfrg_table
    for mse in mses:
        mse.sampid = sampid 
        idx = 1 + add_mse(mse,idx,parent_t,srcfrg_t)
    parent_t.flush()
    srcfrg_t.flush()
    return idx

def save_h5(fname,mses):
    h5t = create_h5_file(fname)
    add_mses(mses,h5t.parent_table,h5t.srcfrg_table)
    return h5t

class MSE(tbls.IsDescription):
    idx = tbls.Int64Col()
    sampid = tbls.StringCol(itemsize=64)
    rt = tbls.Float32Col()
    ccs = tbls.Float64Col()
    mz = tbls.Float64Col()
    ppm = tbls.Float64Col()
    z = tbls.Int8Col()
    intensity = tbls.Float64Col()
    n = tbls.UInt8Col() 
    ms2_data = tbls.Float64Col(shape=(MS2_ARR_LEN,2))
    mgf_files = tbls.StringCol(itemsize=MGF_FILE_STR_LEN)
    src_frag_ids = tbls.UInt64Col(shape=SRC_FRG_ARR_LEN)




