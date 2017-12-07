
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import copy
from matplotlib import cm
from matplotlib.colors import rgb2hex
from graphviz import Digraph
from collections import defaultdict
from .bisect_collection import SortedCollection


H = 1.007825
RT_WINDOW = 0.06#(5/60) # 5sec window aka +/- 2.5sec
CCS_PPT = 10 #10ppt aka 1%
MZ_PPM = 5
MS2_PPM = 25
#print the tolerance values used by the module on import
# can be changed by explicitly overriding in code using this module
# eg.
# >>>import mseutils
# >>>mseutils.MZ_PPM = 10
print("Global Tolerances:")
print("#"*20)
print("RT_WINDOW : {}".format(RT_WINDOW))
print("CCS_PPT : {}".format(CCS_PPT))
print("MZ_PPM : {}".format(MZ_PPM))
print("MS2_PPM : {}".format(MS2_PPM))
print("#"*20)

# color mapper for generating ms2 trees
cmapper = cm.ScalarMappable(cmap='viridis')

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
            # eq1 = (other.val_range[0] <= self.val_range[1] <= other.val_range[1])
            # eq2 = (other.val_range[0] <= self.val_range[0] <= other.val_range[1])
            eqboth = (other.val_range[0] <= self.val_range[1]) and (self.val_range[1] <= other.val_range[1])
            # if eq1 or eq2:
            if eqboth:
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
    
    __slots__ = ['mh', 'z','mz', 'ppm', 'error', 'mz_range','i']
    
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
            try:
                # eq1 = (other.mz_range[0] <= self.mz_range[1] <= other.mz_range[1])
                # eq2 = (other.mz_range[0] <= self.mz_range[0] <= other.mz_range[1])
                eqboth = (other.mz_range[0] <= self.mz_range[1]) and (self.mz_range[1] <= other.mz_range[1])
            except Exception as e:
                print(self.mz_range)
                print(other.mz_range)
                raise e
            # if (eq1 or eq2) and (other.z == self.z):
            if eqboth and (other.z == self.z):
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

    ToDo:
        * use a SortedCollection instead of self.d for faster search and insertion times 
    '''

    def __init__(self,parent_mz,i,iterable=None,filt_thresh=0.01):
        '''
        additional MS2 specific attributes in __init__
        Args:
            parent_mz (MZ) : the mz of the parent ion
            i (float) : intensity of parent ion
            itearable (iterable) : a iterable of 2 tuples in th form [(key,val)]
            filt_thresh (float) : intensity ratio for filtering ms2 ions, see filter()
        '''
        super().__init__(iterable)
        self.parent_mz = parent_mz
        self.filt_thresh = filt_thresh
        self.i = i

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
        dig = Digraph(node_attr={'penwidth':'3'},edge_attr={'penwidth':'2.5'},**kwargs)
        dig.node(name='ParentIon',label=str(self.parent_mz))
        
        items = list(self.items())
        items.sort(key=lambda x:x[0].mz)
        mzs = [x[0] for x in items]
        intensities = [x[1] for x in items]
        colors = cmapper.to_rgba(intensities)

        for i,mzctup in enumerate(zip(mzs,colors)):
            mz,colorv = mzctup
            dig.node(name='DaughterIon_{}'.format(i),label=str(mz),
                color=rgb2hex(colorv[:3]))
            
            dig.edge('ParentIon','DaughterIon_{}'.format(i))

        return dig

    def _repr_svg_(self):
        '''hook for displaying graph in jupyter notebook'''
        dot = self.make_graph()
        binout = dot.pipe('svg')
        return binout.decode('utf-8')

class MSND(MZD):
    '''
    NEEDS (LOTS OF) WORK
    DO NOT USE YET!!!
    Sub-class of MZD Specifically for holding msN data. 
    A few extra methods are defined. 
    '''

    def __init__(self,mz,i,iterable=None,filt_thresh=0.01):
        '''
        additional MS2 specific attributes in __init__
        Args:
            parent_mz (MZ) : the mz of the parent ion
            itearable : a iterable of 2 tuples in th form [(key,val)]
            filt_thresh : intensity ratio for filtering ms2 ions, see filter()
        '''
        super().__init__(iterable)
        self.mz = mz
        self.i = i 
        self.filt_thresh = filt_thresh


    def __contains__(self,other):
        if isinstance(other,self.__class__):
            other = other.mz
        for key in self.d.keys():
            if key == other:
                return True
        return False

    def __getitem__(self,key):
        if isinstance(key,self.__class__):
            key = key.mz
        for skey in self.d.keys():
            if key == skey:
                return self.d[skey]
                
        raise KeyError(str(key))
    
    def _merge(self,msd1,msd2):
        mz1,mz2 = msd1.mz,msd2.mz
        ni = (msd1.i + msd2.i) /2
        nmz = (mz1+mz2) / 2
        new_msd = MSND(mz=nmz,i=ni)
        for k1,v1 in msd1.items():
            if k1 in msd2:
                v2 = msd2[k1]
                nv = (v1 + v2) / 2
                nk = (k1 + msd2[k1]) / 2
                new_msd[nk] = nv
            else:
                new_msd[k1] = v1
        for k2,v2 in msd2.items():
            if k2 not in new_msd:
                new_msd[k2] = v2
        return new_msd



    def __setitem__(self,key,val):
        '''
        looks to see if there is a matching mz within ppm error,
        if one is found then the mz's and intensities are averaged (mean)
        '''
                         
        for skey in self.d.keys():
            if key == skey:
                if isinstance(val,self.__class__): # msN spec
                    nkey = (key+skey) / 2 
                    sval = self.d[skey]
                    if isinstance(sval,self.__class__): #adding to a msN spec
                        val = self._merge(val,sval)
                        nkey.i = val.i
                    else:
                        nkey.i = self.d[skey] #set the intensity to be the ms2 int

                    del self.d[skey]
                    self.d[nkey] = val
                    return

                else: #ms2 spec
                    nkey = (key+skey) / 2 
                    nval = (self.d[skey] + val) / 2
                    del self.d[skey]
                    self.d[nkey] = nval
                    return
        
        if isinstance(val,self.__class__):
            key.i = val.i
        self.d[key] = val #if not found set the val

    def filter(self):
        '''
        function that takes all the current ms2 data and filters ions
        based on intensity. Ions that have intensity/total < the `filt_thresh`
        are thrown out. 
        '''
        total = sum([v if not isinstance(v,self.__class__) else v.mz.i for v in self.values()])
        todel = []
        for key,val in self.items():
            if val / total < self.filt_thresh:
                todel.append(key)
        for key in todel: del self.d[key]

    def __get_layers(self,msnd,layers=[[]],i=0):
        'recursively decend the msn tree and return a list of layers'
        pmz = self.mz
        for k,v in msnd.items():
            if isinstance(v,self.__class__):
                layers[i].append((pmz,k,k.i))
                layers.append([])
                return self._get_layers(v,layers,i+1)
            else:
                layers[i].append((pmz,k,v))
        return layers

    def _get_layers(self,msnd):
        'iteratively decend the ms3 tree and return a list of layers'
        layers = [[],[]]
        pmz = msnd.mz
        for k,v in msnd.items():
            if isinstance(v,self.__class__):
                layers[0].append((pmz,k,k.i))
                pmz2 = v.mz
                for k3,v3 in v.items():
                    layers[1].append((pmz2,k3,v3))
            else:
                layers[0].append((pmz,k,v))
        return layers

    def _add_layer_graph(self,layer,dig):
        'mutator func to add a layer to a digraph'
        intensities = [x[2] for x in layer]
        print(intensities)
        colors = cmapper.to_rgba(intensities)
        dig.node(name=str(layer[0][0]),label=str(layer[0][0]))
        for l,colorv in zip(layer,colors):
            dig.node(name=str(l[1]),label=str(l[1]),
                color=rgb2hex(colorv[:3])) 

            dig.edge(str(l[0]),str(l[1]))

    def make_graph(self,**kwargs):
        '''
        Makes a graphviz graph of the ion tree
        Args:
            **kwargs : keyword args to get passed to `graphviz.Digraph`
        Returns:
            dig : a `graphviz.Digraph` of the ion tree
        '''
        dig = Digraph(node_attr={'penwidth':'3'},edge_attr={'penwidth':'2.5'},**kwargs)
        for layer in self._get_layers(self):
            if layer:
                self._add_layer_graph(layer,dig)
        return dig

    def _repr_svg_(self):
        '''hook for displaying graph in jupyter notebook'''
        dot = self.make_graph()
        binout = dot.pipe('svg')
        return binout.decode('utf-8')

    def __eq__(self,other):
        '''tolerate other ms2d instances for msn style data'''
        raise NotImplemented("can't compare MSND's yet")

class MseSpec(object):
    '''
    This is a class that holds a molecule spectrum from lcms data.
    This could be whole molecule of interest or a source fragment or noise. 
    Contains features like: rt, CCS, and ms2 spectrum as a MZD. 

    '''

    def __init__(self,mz,rt,ccs,ms2_data,i,mgf_files=[],n=0,src_frag=[],rt_window=RT_WINDOW,ccs_ppt=CCS_PPT):
        '''
        Args:
            mz (MZ) : the mz of the parent ion
            rt (RT) : the retention time of the parent ion
            css (CCS) : the ccs of the parent ion
            ms2_data (MS2D) : a MS2D of the ms2 data for the parent
            i (float) : the intensity of the parent ion
            n (int) : number of specs that make up the MseSpec obj
            src_frag (list) : a list of other MseSpecs that are source fragments of this spec 
            mgf_files (list) : a list of the mgf files that gave rise to the data
            rt_window (float) : rt window in minutes
            ccs_window : ccs error in parts per thousand ppt
        ToDo:
            * a lot of possible functionality..
            * add __sub__ methods for blank subtraction?
        '''

        self.mz = mz
        self.rt = RT(rt,window=rt_window)
        self.ccs = CCS(ccs,ppt=ccs_ppt)
        self.ms2_data = ms2_data
        self.mgf_files = mgf_files #needs work
        self.n = n
        self.i = i
        self.src_frags = []

    @classmethod
    def from_dict(cls,indict):
        '''instantiate from a dict, if you are into that kind of thing'''
        return cls(**indict)

    @classmethod
    def from_mgf_dict(cls,mgfd,mgfname):
        '''
        instantiate from a mgf dict generated from pyteomics.mgf.read()
        Args:
            mgfd (dict) : dict from the mgf spec 
            mgfname (str) : the mgf filename
        '''
        def _parse_title(titlestr):
            '''parse the title string from a mgf params dict'''
            tlist = titlestr.split()
            validxs = list(range(1,len(tlist)+1,2))
            vals = [tlist[vi] for vi in validxs]
            keyidxs = list(range(0,len(tlist),2))
            keys = [tlist[ki].rstrip(':') for ki in keyidxs]
            return {k:v for k,v in zip(keys,vals)}


        parent_mass,parent_i = mgfd['params']['pepmass']
        i = float(parent_i)
        tdict = _parse_title(mgfd['params']['title'])
        mz = MZ(
            mz=float(parent_mass),
            z=int(mgfd['params']['charge'][0]))
        
        rt = float(tdict['RetTime'])
        ccs = float(tdict['CCS'])
        mzarr = mgfd['m/z array']
        intarr = mgfd['intensity array']
        MZarr = [MZ(x,ppm=MS2_PPM) for x in mzarr]
        ms2_data = MS2D(mz,i,list(zip(MZarr,intarr)))
        mgf_files = [mgfname]

        return cls(mz=mz, rt=rt, ccs=ccs, i=i,
            ms2_data=ms2_data, mgf_files=mgf_files)


    @property
    def ms2vect(self):
        '''fooling around with a ms2 binned vect'''
        v = np.zeros(2000)
        ms2arr = np.floor(np.array([(ms2mz.mz,i) for ms2mz,i in self.ms2_data.items()]))
        if ms2arr.shape[0] > 0:
            mzvs, inties = ms2arr.T
            relinties = inties / np.sum(inties)
            for mz,i in zip(mzvs,relinties):
                v[int(mz)] += i
            
        return v

    @property
    def n_src_frags(self):
        return len(self.src_frags)
    
    def add_src_frag(self,other):
        self.src_frags.append(other)

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
        if not isinstance(other,self.__class__):
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
    
    # def __repr__(self):
    #     return "{classn}(mz={mz},rt={rt},ccs={ccs},{ms2_data})".format(
    #         classn = self.__class__.__name__,
    #         mz = self.mz,
    #         rt = self.rt,
    #         ccs = self.ccs,
    #         ms2_data = self.ms2_data)
    def __repr__(self):
        return "{classn}(mz={mz},rt={rt},ccs={ccs})".format(
            classn = self.__class__.__name__,
            mz = self.mz,
            rt = self.rt,
            ccs = self.ccs)


    ########
    # Hash #  
    ########

    def __hash__(self):
        return hash(self.__repr__())


def src_frags(mol_specs):
    '''
    Takes a list of combined MseSpecs and looks for source fragments
    This is done by looking in a given rt window for MseSpec's whose parent
    ion is in the ms2_data of another spec.
    '''
    sms = SortedCollection(mol_specs,key=lambda x:x.rt.val)
    src_frg_idxs = defaultdict(list)
    for i,ms1 in enumerate(tqdm(sms)):
        rt_chunk = sms.find_between(*ms1.rt.val_range)
        if rt_chunk:
            for ms2 in rt_chunk:
                if ms2.mz in ms1.ms2_data:
                    src_frg_idxs[i].append(ms2)
                else:
                    src_frg_idxs[i] = []
        else:
            # tqdm.write("No rt elements between {} and {}".format(*ms1.rt.val_range)) # debug line
            src_frg_idxs[i] = []
    print("combining srg frags...")
    combined = []
    for idx,frgs in src_frg_idxs.items():
        parent = sms[idx]
        for frg in frgs:
            if frg is not parent:
                parent.add_src_frag(frg)
        combined.append(parent)
    return combined
