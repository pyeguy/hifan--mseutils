from pyteomics import mgf
import os
from .utils import *

def load_mgf(mgfname):
    '''load a mgf file and output the MseSpecs'''
    mgf_reader = mgf.read(mgfname)
    filename = os.path.basename(mgfname).split('.')[0]
    # commented out in favor of tqdm for loop
    # mses = [MseSpec.from_mgf_dict(spec,filename) for spec in mgf_reader]
    mses = []
    # for spec in tqdm(mgf_reader,desc='loading mgf'):
    for spec in mgf_reader: # no pbar
        mses.append(MseSpec.from_mgf_dict(spec,filename))        
    return mses


def load_frag_csv(csv_file,mz_kwargs={},molspec_kwargs={}):
    '''
    Loads a csv file of that includes fragments and returns a list of MseSpec
    '''
    # print("reading csv file and grouping...")
    df = pd.read_csv(csv_file)
    # df = df[df.Ar1 > 0.33]
    gb = df.groupby(("RetTime","CCS","PrecMz","PrecZ","MgfFileName",'PrecIntensity'))
    mss = []
    # print('creating MseSpecs')
    # for gtup,gdf in tqdm(gb,"processing file"):
    for gtup,gdf in gb: # no pbar
        rt,ccs,mz,z,mgf_fname,preci = gtup
        mzo = MZ(mz=mz,z=z,**mz_kwargs)

        ms2vals = gdf[['ProdMz','ProdIntensity']].values
        if ms2vals.any():
            ms2 = MS2D(mz,preci,[(MZ(mz,ppm=MS2_PPM),i) for mz,i in ms2vals])
            ms2.filter()
        else:
            ms2 = MS2D(mz,preci)
        ms = MseSpec(
            mz=mzo,
            rt=rt,
            ccs=ccs,
            ms2_data=ms2,
            mgf_files=[mgf_fname],
            i = preci, 
            **molspec_kwargs)
        mss.append(ms)
    return mss

def load_rep_csv(csv_file):
    '''
    Loads a csv file of replicate MS1 masses and returns a list 
    of MseSpec's 
    '''
    df = pd.read_csv(csv_file)
    df['z'] = df.Precursor.apply(lambda x:x[1])
    mss = []
    for idx,row in df.iterrows():
        mz = MZ(row["PrecMz"],row["z"])
        ccs = row["CCS"]
        rt = row["RetTime"]
        i = row["PrecIntensity"]
        ms2_data = MS2D(mz,i)
        ms = MseSpec(
            mz=mz,
            ccs=ccs,
            rt=rt,
            ms2_data=ms2_data,
            i=i)
        mss.append(ms)
    return mss

def load_rep_and_frags_csv(rep_csv,frag_csv_file,mz_kwargs={},msespec_kwargs={}):
    '''
    Loads a MS1 rep file and the MS2 rep_frag file and then combines all the frags
    inso the MS1 masses where they exist.
    Args:
        rep_csv (str) : csv file of the replicated MS1's
        frag_csv (str) : csv file of the rep fragments
        mz_kwargs (dict) : a dict of any kwargs to pass to the MZ instantiation
        msespec_kwargs (dict) : a dict of any kwargs to pass to the MseSpec instantiation
    Returns:
        cmss (list) : combined replicate MseSpecs 
    '''
    pmss = load_rep_csv(rep_csv)
    mss = load_frag_csv(frag_csv_file)
    smss = SortedCollection(mss,key= lambda x:x.rt.val)
    cmss = []

    # for search_ms in tqdm(pmss,desc="combining MseSpec's"):
    for search_ms in pmss: # no pbar
        search_ms = copy(search_ms) # added copy here... necissary?
        rt_chunk = smss.find_between(*search_ms.rt.val_range)
        if rt_chunk:
            for i,ms in enumerate(rt_chunk):
                try:
                    if search_ms == ms:
                        search_ms += ms

                except Exception as e:
                    print(search_ms.__repr__())
                    print(ms.__repr__())
                    raise e
        if len(search_ms.mgf_files) >= 2:
            cmss.append(search_ms)

    # comb = len(cmss)
    # print("{} combined spec ({:.2f}%)".format(comb,(len(mss)-comb)/len(mss) *100))
    return cmss

def load_h5(fname,group_name='msedata',table_name='mse_specs'):
    '''dont forget to add in source frags'''
    pass