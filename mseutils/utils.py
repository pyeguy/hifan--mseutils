
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import copy

from graphviz import Digraph


H = 1.007825
RT_WINDOW = 0.6#(5/60) # 5sec window aka +/- 2.5sec
CCS_PPT = 10 #10ppt aka 1%
MZ_PPM = 5
MS2_PPM = 25

class Chunker(object):
    '''
    this is a work in progress for making a chunking class.
    ToDo: 
        * make it a true generator class 
            * think of a way to guess length 
    '''
    def __init__(self,tochunk,chunk_size=250):
        self.tochunk = tochunk
        self.chunk_size = chunk_size
        self.cnt = 0
    
    def __len__(self):
        if len(self.tochunk) % self.chunk_size == 0:
            return len(self.tochunk) // self.chunk_size
        else:
            return (len(self.tochunk) // self.chunk_size) + 1
    
    def __iter__(self):
        return self # see __next__ for iteration

    def __next__(self):
        self.cnt +=1
        if self.cnt > len(self):
            raise StopIteration
        i = (self.cnt  - 1) * self.chunk_size
        try:
            return self.tochunk[i:i+self.chunk_size]
        except IndexError:
            return self.tochunk[i:]

class FuzzyCompare(object):
    '''
    Base class for comparing things with error ranges.
    '''
    __slots__ = ('val','val_range','window','error','error_func')
    
    def __init__(self,val,window=0,error_func=None):
        '''
        Args:
            val (number) : the value of the object
            window (number) : error window aka val +/- window/2
            error_func (function) : a function for computing error, used for ppx style error
        '''
        self.val = val
        self.error_func = error_func
        if self.error_func is not None:
            self.error =  error_func(self.val)
            self.window = 2 * self.error
        else:
            self.window = window
            self.error = window / 2
        self.val_range = (self.val - (self.window /2), self.val+ (self.window/2))

    ################
    # Comparisons  #
    ################
    
    def __lt__(self, other):
        if isinstance(other,self.__class__):
            other = other.val_range[0]

        if self.val_range[1] < other:
            return True
        else:
            return False

    def __le__(self, other):
        if self.__lt__(other) or self.__eq__(other):
            return True
        else:
            return False

    def __eq__(self, other):
        if isinstance(other,self.__class__):
            eq1 = (other.val_range[0] <= self.val_range[1] <= other.val_range[1])
            eq2 = (other.val_range[0] <= self.val_range[0] <= other.val_range[1])
            if eq1 or eq2:
                return True
            else:
                return False
        else:
            if self.val_range[0] <= other <= self.val_range[1]:
                return True
            else:
                return False
                
    def __ne__(self, other):
        if not self.__eq__(other):
            return True
        else:
            return False

    def __gt__(self, other):
        if isinstance(other,self.__class__):
            other = other.val_range[1]
        
        if self.val_range[0] > other:
            return True
        else:
            return False

    def __ge__(self, other):
        if self.gt(other) or self.__eq__(other):
            return True
        else:
            return False

    ########
    # Math #
    ########

    def __add__(self, other):
        if isinstance(other,self.__class__):
            nval = self.val + other.val
            if self.error_func and other.error_func:
                # get earch error for a random number (100)
                serror = self.error_func(100)
                oerror = other.error_func(100)
                # pick the function that gives the largest error
                nerror_func = self.error_func if serror >= oerror else other.error_func
                return self.__class__(nval,error_func=nerror_func)
            else:
                nwindow = self.window if self.window >= other.window else other.window
                return self.__class__(nval,window=nwindow)
        else:
            raise TypeError(type(other))
    
    def __radd__(self,other):
        return self.__add__(other)

    def __truediv__(self, other):
        if isinstance(other,self.__class__):
            nval = self.val + other.val
            if self.error_func and other.error_func:
                # get earch error for a random number (100)
                serror = self.error_func(100)
                oerror = other.error_func(100)
                # pick the function that gives the largest error
                nerror_fucnc = self.error_func if serror >= oerror else other.error_func
                return self.__class__(nval,error_func=nerror_fucnc)
            else:
                nwindow = self.window if self.window >= other.window else other.window
                return self.__class__(nval,window=nwindow)
        else:
            nval = self.val / other
            if self.error_func is not None:
                return self.__class__(nval,error_func=self.error_func)
            else:
                return self.__class__(nval,window=self.window)



    def __str__(self):
        return "{val:.4f} +/- {error:.2g}".format(val=self.val, error=self.error)

class RT(FuzzyCompare):
    '''Retention time pass through class of FuzzyCompare'''
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

class CCS(FuzzyCompare):
    '''FuzzyCompare subclass for ccs values'''
    def __init__(self,val,ppt=0,error_func=None):
        '''
        Args:
            val (flaot) : ccs value
            ppt (float) : error in parts per thousand, used if error_func isn't specified
            error_func (func): an error function for ccs error (should be probably be monotonic w/ val)
        '''
        if error_func is None:
            error_func = lambda x:x * (ppt/1e3)

        super().__init__(val,error_func=error_func)
    
    def __str__(self):
        return "{val:.4f} +/- {error:.2g}".format(val=self.val, error=self.error)

class MZ(object):
    """Numerical class to hold m/z values and do rich comparisons 
    and (some) math using the ppm error ranges.
    
    ToDo:
        * refactor to be a subclass of FuzzyCompare (i wrote this first..)

    Notes: 
    All math is currently done on the `mz` value, maybe it
    should be on the `m` value. I don't really know when you would
    use math on these objects but it feels like you might...


    """
    
    __slots__ = ['mh', 'z','mz', 'ppm', 'error', 'mz_range']
    
    def __init__(self,mz,z=1,ppm=MZ_PPM):
        '''
        Can be instantiated using just m if z = 1 and ppm = 5.
        
        Args:
            mz : precision mz
            z : charge (defaults to 1 making mz = m)
            ppm : ppm error/precision of instrument 

        '''
        self.mz = mz
        self.z = int(z)
        self.ppm = ppm
        self.mh = (self.mz / self.z) - ((self.z-1) * H) #assumes only H adducts
        self.error = self.mz * (self.ppm / 10e6)
        self.mz_range = (
            (self.mz - self.error), 
            (self.mz + self.error))

    ################
    # Comparisons  #
    ################
    
    def __lt__(self, other):
        if isinstance(other,MZ):
            other = other.mz_range[0]

        if self.mz_range[1] < other:
            return True
        else:
            return False

    def __le__(self, other):
        if self.__lt__(other) or self.__eq__(other):
            return True
        else:
            return False

    def __eq__(self, other):
        if isinstance(other,MZ):
            eq1 = (other.mz_range[0] <= self.mz_range[1] <= other.mz_range[1])
            eq2 = (other.mz_range[0] <= self.mz_range[0] <= other.mz_range[1])
            if (eq1 or eq2) and (other.z == self.z):
                return True
            else:
                return False
        else:
            if self.mz_range[0] <= other <= self.mz_range[1]:
                return True
            else:
                return False
                
    def __ne__(self, other):
        if not self.__eq__(other):
            return True
        else:
            return False

    def __gt__(self, other):
        if isinstance(other,MZ):
            other = other.mz_range[1]
        
        if self.mz_range[0] > other:
            return True
        else:
            return False

    def __ge__(self, other):
        if self.gt(other) or self.__eq__(other):
            return True
        else:
            return False

    ########
    # Math #
    ########

    def __add__(self, other):
        if isinstance(other,MZ):
            oppm = other.ppm
            other = other.mz
            if oppm > self.ppm:
                ppm = oppm
            else:
                ppm = self.ppm
        else:
            ppm = self.ppm

        return (MZ(self.mz + other,ppm=ppm))
    
    def __radd__(self,other):
        return self.__add__(other)


    def __sub__(self, other):
        if isinstance(other,MZ):
            oppm = other.ppm
            other = other.mz
            if oppm > self.ppm:
                ppm = oppm
            else:
                ppm = self.ppm
        else:
            ppm = self.ppm

        return (MZ(self.mz - other,ppm=ppm))

    def __rsub__(self,other):
        return self.__sub__(other)


    def __mul__(self, other):
        if isinstance(other,MZ):
            raise NotImplemented
        else:
            return MZ(self.mz * other,ppm=self.ppm)
    
    def __rmul__(self,other):
        return self.__mul__(other)


    def __matmul__(self, other):
        raise NotImplemented

    def __truediv__(self, other):
        if isinstance(other,MZ):
            raise NotImplemented
        else:
            return MZ(self.mz / other,ppm=self.ppm)

    def __floordiv__(self, other):
        if isinstance(other,MZ):
            raise NotImplemented
        else:
            return MZ(self.mz // other,ppm=self.ppm)

    def __mod__(self, other):
        if isinstance(other,MZ):
            raise NotImplemented
        else:
            return MZ(self.mz % other,ppm=self.ppm)

    def __divmod__(self, other):
        raise NotImplemented

    def __pow__(self, other):
        if isinstance(other,MZ):
            raise NotImplemented
        else:
            return MZ(self.mz ** other,ppm=self.ppm)

    ##########################
    # String Representations #
    ##########################

    def __str__(self):
        return "{mz:.4f} +/- {error:.2g}".format(mz=self.mz, error=self.error)

    def __repr__(self):
        return "{classn}(mz={mz},z={z},ppm={ppm})".format(
            classn = self.__class__.__name__,
            mz = self.mz,
            z = self.z,
            ppm = self.ppm)

    ########
    # Hash #  
    ########

    def __hash__(self):
        return hash(self.__repr__())

class MZD(object):
    '''
    dict like object for storing {MZ:intensity} data.
    overrides __contains__ in order to do MZ __eq__ method
    Note:
    MZ objects are mutable, weird things can happen making dicts from 
    mutable ojbects... might re-implement in Cython to solve this problem. 
    see below for example:
        https://stackoverflow.com/questions/4828080/how-to-make-an-immutable-object-in-python

    '''
    def __init__(self,iterable=None):
        '''
        Args:
            iterable : an iterable of 2 tuples of the form [(key,val)]
        '''
        self.d = dict()
        if iterable:
            for key,val in iterable:
                self[key]=val

    def __contains__(self,other):
        for key in self.d.keys():
            if key == other:
                return True
        return False

    def __getitem__(self,key):
        for skey in self.d.keys():
            if key == skey:
                return self.d[skey]

        raise KeyError(str(key))
    
    def __setitem__(self,key,val):
        '''
        looks to see if there is a matching mz within ppm error,
        if one is found then the mz's and intensities are averaged (mean)
        '''
        for skey in self.d.keys():
            if key == skey:
                nkey = (key+skey) / 2 
                nval = (self.d[skey] + val) / 2
                del self.d[skey]
                self.d[nkey] = nval
                return
        self.d[key] = val

    def update(self,other_dict):
        '''equivalent functionality of dict.update but using custom __setitem__'''
        for key,val in other_dict.items():
            self[key] = val
    
    def items(self):
        for key,val in self.d.items():
            yield key,val
    def keys(self):
        return self.d.keys()
    
    def values(self):
        return self.d.values()
    
    def __iter__(self):
        return self.d.__iter__()

    def __len__(self):
        return len(self.d)

    def __repr__(self):
        return self.d.__repr__()

    def __str__(self):
        return str(self.d)

class MS2D(MZD):
    '''
    Sub-class of MZD Specifically for holding ms2 data. 
    A few extra methods are defined. 
    '''

    def __init__(self,parent_mz,iterable=None,filt_thresh=0.01):
        '''
        additional MS2 specific attributes in __init__
        Args:
            parent_mz (MZ) : the mz of the parent ion
            itearable : a iterable of 2 tuples in th form [(key,val)]
            filt_thresh : intensity ratio for filtering ms2 ions, see filter()
        '''
        super().__init__(iterable)
        self.parent_mz = parent_mz
        self.filt_thresh = filt_thresh

    def filter(self):
        '''
        function that takes all the current ms2 data and filters ions
        based on intensity. Ions that have intensity/total < the `filt_thresh`
        are thrown out. 
        '''
        total = sum(self.values())
        todel = []
        for key,val in self.items():
            if val / total < self.filt_thresh:
                todel.append(key)
        for key in todel: del self.d[key]


    def make_graph(self,**kwargs):
        '''
        Makes a graphviz graph of the ion tree
        Args:
            **kwargs : keyword args to get passed to `graphviz.Digraph`
        Returns:
            dig : a `graphviz.Digraph` of the ion tree
        '''
        dig = Digraph(**kwargs)
        dig.node(name='ParentIon',label=str(self.parent_mz))
        for i,mz in enumerate(self.keys()):
            dig.node(name='DaughterIon_{}'.format(i),label=str(mz))
            dig.edge('ParentIon','DaughterIon_{}'.format(i))
        return dig

    def _repr_svg_(self):
        '''hook for displaying graph in jupyter notebook'''
        dot = self.make_graph()
        binout = dot.pipe('svg')
        return binout.decode('utf-8')

class MseSpec(object):
    '''
    This is a class that holds a molecule spectrum from lcms data.
    This could be whole molecule of interest or a source fragment or noise. 
    Contains features like: rt, CCS, and ms2 spectrum as a MZD. 

    '''

    def __init__(self,mz,rt,ccs,ms2_data,i,mgf_files=[],n=0,rt_window=RT_WINDOW,ccs_ppt=CCS_PPT):
        '''
        Args:
            mz (MZ) : the mz of the parent ion
            rt (RT) : the retention time of the parent ion
            css (CCS) : the ccs of the parent ion
            ms2_data (MS2D) : a MS2D of the ms2 data for the parent
            i (float) : the intensity of the parent ion
            n (int) : number of specs that make up the MseSpec obj
            mgf_files (list) : a list of the mgf files that gave rise to the data
            rt_window (float) : rt window in minutes
            ccs_window : ccs error in parts per thousand ppt
        ToDo:
            * a lot...
            * add __sub__ methods for blank subtraction?
        '''

        self.mz = mz
        self.rt = RT(rt,window=rt_window)
        self.ccs = CCS(ccs,ppt=ccs_ppt)
        self.ms2_data = ms2_data
        self.mgf_files = mgf_files #needs work
        self.n = n
        self.i = i


    @classmethod
    def from_dict(cls,indict):
        '''instantiate from a dict, if you are into that kind of thing'''
        return cls(**indict)

    @property
    def ms2vect(self):
        '''fooling around with a ms2 binned vect'''
        v = np.zeros(2000)
        ms2arr = np.floor(np.array([(ms2mz.mz,i) for ms2mz,i in self.ms2_data.items()]))
        for mz,i in ms2arr:
            v[int(mz)] += i
        return v

    ################
    # Comparisons  #
    ################
   
    def __eq__(self, other):
        '''
        checks mz, rt and ccs equality right now, ms2 tbd
        '''
        if isinstance(other,self.__class__):
            rteq = self.rt == other.rt
            ccseq = self.ccs == other.ccs
            mzeq = self.mz == other.mz

            if rteq and ccseq and mzeq:
                return True
            else:
                return False
        else:
            raise TypeError("Can only compare to other MseSpec instances")
                
    def __ne__(self, other):
        if not self.__eq__(other):
            return True
        else:
            return False
    
    def __add__(self,other):
        if not isinstance(other,self):
            raise TypeError(type(other))
        nmz = (self.mz + other.mz) / 2
        nrt = (self.rt + other.rt) / 2
        nccs = (self.ccs + other.ccs) / 2 
        nms2_data = copy(self.ms2_data).update(other.ms2_data)
        nn = self.n + 1
        i = (self.i + other.i) / 2 
        mgf_files = self.mgf_files + other.mgf_files
        return MseSpec(nmz,nrt,nccs,nms2_data,n=nn,i=i,mgf_files=mgf_files)

    def __radd__(self,other):
        return self.__add__(self,other)

    def __iadd__(self,other):
        if not isinstance(other,self.__class__):
            raise TypeError(type(other))
        self.mz = (self.mz + other.mz) / 2
        self.rt = (self.rt + other.rt) / 2
        self.ccs = (self.ccs + other.ccs) / 2 
        self.ms2_data.update(other.ms2_data)
        self.i = (self.i + other.i) / 2
        self.n += 1
        self.mgf_files.extend(other.mgf_files)
        return self


    ##########################
    # String Representations #
    ##########################

    def __str__(self):
        outstr = "mz : {mz} @ rt: {rt} & ccs : {ccs}".format(
            mz=self.mz,
            rt=self.rt,
            ccs=self.ccs)
        return outstr
    
    def __repr__(self):
        return "{classn}(mz={mz},rt={rt},ccs={ccs},{ms2_data})".format(
            classn = self.__class__.__name__,
            mz = self.mz,
            rt = self.rt,
            ccs = self.ccs,
            ms2_data = self.ms2_data)

    ########
    # Hash #  
    ########

    def __hash__(self):
        return hash(self.__repr__())




def load_csv(csv_file,mz_kwargs={},molspec_kwargs={}):
    '''
    Loads a csv file of frags and returns a list of MseSpec
    '''
    df = pd.read_csv(csv_file)
    # df = df[df.Ar1 > 0.33]
    gb = df.groupby(("RetTime","CCS","PrecMz","PrecZ","MgfFileName",'PrecIntensity'))
    mss = []
    for gtup,gdf in tqdm(gb,"loading file"):
        rt,ccs,mz,z,mgf_fname,preci = gtup
        mzo = MZ(mz=mz,z=z,**mz_kwargs)

        ms2vals = gdf[['ProdMz','ProdIntensity']].values
        if ms2vals.any():
            ms2 = MS2D(mz,[(MZ(mz,ppm=MS2_PPM),i) for mz,i in ms2vals])
            ms2.filter()
        else:
            ms2 = MS2D(mz)
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
        ms2_data = MS2D(mz)
        i = row["PrecIntensity"]
        ms = MseSpec(
            mz=mz,
            ccs=ccs,
            rt=rt,
            ms2_data=ms2_data,
            i=i)
        mss.append(ms)
    return mss



def load_rep_frags_csv(rep_csv,frag_csv_file,mz_kwargs={},msespec_kwargs={}):
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
    mss = load_csv(frag_csv_file)
    orig_len = len(mss)
    worst_case = (orig_len**2//2)-orig_len
    pbar = tqdm(total=worst_case,desc="combining MseSpec's")
    cmss = []

    while len(pmss) > 0:
        search_ms = pmss.pop(0)
        if len(mss) > 0:
            todel = []
            for i,ms in enumerate(mss):
                try:
                    if search_ms == ms:
                        search_ms += ms
                        todel.append(i)
                except Exception as e:
                    print(search_ms.__repr__())
                    print(ms.__repr__())
                    raise e
            todel.sort(reverse=True)
            for i in todel: del mss[i]
        if len(search_ms.mgf_files) >= 2:
            cmss.append(search_ms)
        pbar.update(orig_len - len(mss))
    pbar.total = pbar.n #set the prog bar to 100%
    pbar.close()
    comb = len(cmss)
    print("{} combined spec ({:.2f}%)".format(comb,(orig_len-comb)/orig_len *100))
    return cmss



## old func.. should go but keeping it around for sentimental reasons
# def load_replicate_fragments_csv(csv_file,mz_kwargs={},molspec_kwargs={}):
#     mss = load_csv(csv_file, mz_kwargs, molspec_kwargs)
#     orig_len = len(mss)
#     worst_case = (orig_len**2//2)-orig_len
#     pbar = tqdm(total=worst_case,desc="combining MseSpec's")
#     cmss = []
#     while len(mss) > 0:
#         search_ms = mss.pop(0)
#         if len(mss) > 0:
#             todel = []
#             for i,ms in enumerate(mss):
#                 try:
#                     if search_ms == ms:
#                         search_ms += ms
#                         todel.append(i)
#                 except Exception as e:
#                     print(search_ms.__repr__())
#                     print(ms.__repr__())
#                     raise e
#             todel.sort(reverse=True)
#             for i in todel: del mss[i]
#         if len(search_ms.mgf_files) >= 2:
#             cmss.append(search_ms)
#         pbar.update(orig_len - len(mss))
#     pbar.total = pbar.n #set the prog bar to 100%
#     pbar.close()
#     comb = len(cmss)
#     print("{} combined spec ({:.2f}%)".format(comb,(orig_len-comb)/orig_len *100))
#     return cmss


## just a way to sanity check the non pw version
# import networkx as nx
# def load_replicate_csv_pw(csv_file,mz_kwargs={},molspec_kwargs={}):
#     mss = load_csv(csv_file, mz_kwargs, molspec_kwargs)
#     orig_len = len(mss)
#     worst_case = (orig_len**2//2)-orig_len
#     pbar = tqdm(total=worst_case,desc="combining MseSpec's")
#     cmss_idx = []
#     cmss = []
#     for i,ms1 in enumerate(mss):
#         for j,ms2 in enumerate(mss[i+1:]):
#             if ms1 == ms2:
#                 cmss_idx.append((i,i+j+1))

#             pbar.update()
#     pbar.close()
#     seen = set([idx for tup in cmss_idx for idx in tup])
#     g = nx.Graph()
#     g.add_edges_from(cmss_idx)
#     graphs = graphs = list(nx.connected_component_subgraphs(g)) 
#     print("Adding Graphs")
#     for graph in graphs:
#         nodes = list(graph.nodes())
#         to_grow = mss[nodes.pop(0)]
#         for node in nodes:
#             to_grow += mss[node]
#         cmss.append(to_grow)
#     unique = [ms for i,ms in enumerate(mss) if i not in seen]
#     cmss.extend(unique)

#     comb = len(cmss)
#     print("{} combined spec ({:.2f}%)".format(comb,(orig_len-comb)/orig_len *100))
#     return cmss