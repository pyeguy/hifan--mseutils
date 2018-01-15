from pyteomics import mgf
import os
from .utils import *
from . import mseh5

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
            mgf_files={mgf_fname},
            i = preci, 
            **molspec_kwargs)
        mss.append(ms)
    return mss


def load_frag_csv_from_prod_seq(csv_file,forceRep=False,mz_kwargs={},molspec_kwargs={}):
    '''
    Loads a csv file of that includes fragments and returns a list of MseSpec
    '''
    # print("reading csv file and grouping...")
    df = pd.read_csv(csv_file)
    # df = df[df.Ar1 > 0.33]
    gb = df.groupby(("MgfFileName","Precursor"))
    combd = defaultdict(list)
    # print('creating MseSpecs')
    # for gtup,gdf in tqdm(gb,"processing file"):
    
    for gtup,gdf in gb: # no pbar
        bad_egg = False
        mgf_fname,prec, = gtup
        kwargs = {"Precursor":prec}
        kwargtups = (
                    ('rt',"RetTime"),
                    ('ccs',"CCS"),
                    ('mz',"PrecMz"),
                    ('z',"PrecZ"),
                    ('i','PrecIntensity'),
                    ('sampid','Sample'))
        for key,valstr in kwargtups:
            valset = set(gdf[valstr])
            if len(valset) !=1 :
                bad_egg = True
                # print(valstr,':',valset)
                # raise ValueError('Too many unique {} values.'.format(valstr))
                # raise ValueError('Too many unique {} values.'.format(valstr))

            val = valset.pop()
            kwargs[key] = val

        if bad_egg:
            continue

        mzo = MZ(mz=kwargs['mz'],z=kwargs['z'],**mz_kwargs)
        del kwargs['mz']
        del kwargs['z']

        ms2vals = gdf[['ProdMz','ProdIntensity']].values
        if ms2vals.any():
            ms2 = MS2D(mzo,kwargs['i'],[(MZ(mz,ppm=MS2_PPM),i) for mz,i in ms2vals])
            ms2.filter()
        else:
            ms2 = MS2D(mzo,kwargs['i'])
        
        mse = MseSpec(
            mz=mzo,
            ms2_data=ms2,
            mgf_files={mgf_fname},
            **kwargs)
        combd[mse.Precursor].append(mse)

    if forceRep:
        todel = [k for k,v in combd.items() if len(v) < 2]
        for k in todel: del combd[k]

    cmss = [_quick_comb(mses) for mses in combd.values()]
    return cmss

def _quick_comb(mses):
    parent = mses.pop(0)
    for mse in mses:
        parent += mse
    return parent


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

def load_rep_and_frags_csv(rep_csv,frag_csv_file,mz_kwargs={},msespec_kwargs={},forceRep=False):
    '''
    Loads a MS1 rep file and the MS2 rep_frag file and then combines all the frags
    inso the MS1 masses where they exist.
    Args:
        rep_csv (str) : csv file of the replicated MS1's
        frag_csv (str) : csv file of the rep fragments
        mz_kwargs (dict) : a dict of any kwargs to pass to the MZ instantiation
        msespec_kwargs (dict) : a dict of any kwargs to pass to the MseSpec instantiation
    Returns:
        combined_mss (list) : combined replicate MseSpecs 
    '''
    parent_mss = load_rep_csv(rep_csv)
    frag_mses = load_frag_csv(frag_csv_file)
    sc_frag_mses = SortedCollection(frag_mses,key= lambda x:x.rt.val)
    combined_mss = []
    
    # for search_ms in tqdm(pmss,desc="combining MseSpec's"):
    for search_ms in parent_mss: # no pbar
        # search_ms = copy(search_ms) # added copy here... necissary?
        rt_chunk = sc_frag_mses.find_between(*search_ms.rt.val_range)
        if rt_chunk:
            for i,ms in enumerate(rt_chunk):
                if search_ms == ms:
                    try:
                        search_ms += ms
                    except Exception as e:
                        print("Exception in Combining following MseSpecs:")
                        print(search_ms.__repr__())
                        print(ms.__repr__())
                        raise e
        if forceRep:
            if len(search_ms.mgf_files) > 1:
                combined_mss.append(search_ms)
        else:
            combined_mss.append(search_ms)
    # comb = len(combined_mss)
    # print("{} combined spec ({:.2f}%)".format(comb,(len(mss)-comb)/len(mss) *100))
    return combined_mss

def load_h5(fname,mode='r', group_name='msedata',parent_tbl='mse_specs', srcfrg_tbl='source_frags',head=None):
    '''dont forget to add in source frags'''
    h5t = mseh5.open_h5_file(fname,mode,group_name,parent_tbl,srcfrg_tbl)
    parent_table = h5t.parent_table
    srcfrg_table = h5t.srcfrg_table
    mses = []
    for i,row in enumerate(tqdm(parent_table.iterrows(),total=len(parent_table))):
        if head:
            if i>head:
                break
        mse = MseSpec.from_h5_row(row)
        if mse.src_frag_ids:
            frg_ids = sorted([x for x in mse.src_frag_ids if x!=0])
            mini,maxi = frg_ids[0],frg_ids[-1]
            frags = srcfrg_table.where("({mini} <= idx) & (idx <= {maxi})".format(mini=mini,maxi=maxi))
            for row in frags:
                mse.add_src_frag(MseSpec.from_h5_row(row))
        mses.append(mse)

    return mses

