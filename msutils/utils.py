
import pandas as pd
from tqdm import tqdm

H = 1.007825

class FuzzyCompare(object):
    __slots__ = ('val','val_range','window','error')
    
    def __init__(self,val,window=0,error_func=None):
        self.val = val
        if error_func is not None:
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
        if isinstance(other,FuzzyCompare):
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
        if isinstance(other,FuzzyCompare):
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
        if isinstance(other,FuzzyCompare):
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
    def __str__(self):
        return "{val:.4f} +/- {error:.2g}".format(val=self.val, error=self.error)

class RT(FuzzyCompare):
    pass

class CCS(FuzzyCompare):
    def __init__(self,val,ppt):
        error_func = lambda x:x * (ppt/1e3)
        super().__init__(val,error_func=error_func)
    
    def __str__(self):
        return "{val:.4f} +/- {error:.2g}".format(val=self.val, error=self.error)



class MZ(object):
    """Numerical class to hold m/z values and do rich comparisons 
    and (some) math using the ppm error ranges.
    
    Note: All math is currently done on the `mz` value, maybe it
    should be on the `m` value. I don't really know when you would
    use math on these objects but it feels like you might...
    """
    
    __slots__ = ['mh', 'z','mz', 'ppm', 'error', 'mz_range']
    
    def __init__(self,mz,z=1,ppm=5):
        '''
        Can be instantiated using just m if z = 1 and ppm = 5.
        
        Args:
            mz : precision mz
            z : charge (defaults to 1 making mz = m)
            ppm : ppm error/precision of instrument (default = 5)

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
            if eq1 or eq2:
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

class MZD(dict):
    '''
    dict like object for storing {MZ:intensity} data.
    overrides __contains__ in order to do MZ __eq__ method
    Note:
    MZ objects are mutable, weird things can happen making dicts from 
    mutable ojbects... might re-implement in Cython to solve this problem. 
    see below for example:
        https://stackoverflow.com/questions/4828080/how-to-make-an-immutable-object-in-python

    '''
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
    def __contains__(self,other):
        for key in self.keys():
            if key == other:
                return True
        return False

class MolSpec(object):
    '''
    This is a class that holds a molecule spectrum from lcms data.
    This could be whole molecule of interest or a source fragment or noise. 
    Contains features like: rt, CCS, and ms2 spectrum as a MZD. 

    '''

    def __init__(self,mz,rt,ccs,ms2_data,rt_window=(5/60),ccs_ppt=5):
        self.mz = mz
        self.rt = RT(rt,window=rt_window)
        self.ccs = CCS(ccs,ppt=ccs_ppt)
        self.ms2_data = ms2_data


    @classmethod
    def from_dict(cls,indict):
        rt = indict['rt']
        css = indict['css']
        ms2_data = indict['ms2_data']
        return cls(rt,css,ms2_data)

    ################
    # Comparisons  #
    ################
    

    def __eq__(self, other):
        '''
        checks mz, rt and ccs equality right now, ms2 tbi
        '''
        if isinstance(other,MolSpec):
            rteq1 = (other.rt_range[0] <= self.rt_range[1] <= other.rt_range[1])
            rteq2 = (other.rt_range[0] <= self.rt_range[0] <= other.rt_range[1])
            rteq = rteq1 or rteq2

            ccseq1 = (other.ccs_range[0] <= self.ccs_range[1] <= other.ccs_range[1])
            ccseq2 = (other.ccs_range[0] <= self.ccs_range[0] <= other.ccs_range[1])
            ccseq = ccseq1 or ccseq2

            mzeq = self.mz == other.mz

            if rteq and ccseq and mzeq:
                return True
            else:
                return False
        else:
            raise TypeError("Can only compare to other MolSpec instances")
                
    def __ne__(self, other):
        if not self.__eq__(other):
            return True
        else:
            return False

    def __str__(self):
        outstr = "mz : {mz} @ rt: {rt} & ccs : {ccs}".format(
            mz=self.mz,
            rt=self.rt,
            ccs=self.ccs)
        return outstr
    
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


def load_csv(csv_file,mz_kwargs={},molspec_kwargs={}):
    df = pd.read_csv(csv_file)
    gb = df.groupby(("RetTime","CCS","PrecMz","PrecZ"))
    mss = []
    for gtup,gdf in tqdm(gb):
        rt,ccs,mz,z = gtup
        mzo = MZ(mz=mz,z=z,**mz_kwargs)
        ms2vals = gdf[['ProdMz','ProdIntensity']].values
        ms2 = MZD([(MZ(mz),i) for mz,i in ms2vals])
        # for idx,row in gdf.iterrows():
        #     ms2mz = MZ(mz=row['ProdMz'])
        #     i = row['ProdIntensity']
        #     ms2[ms2mz]=i
        ms = MolSpec(
            mz=mzo,
            rt=rt,
            ccs=ccs,
            ms2_data=ms2,
            **molspec_kwargs)
        mss.append(ms)
    return mss



