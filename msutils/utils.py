

class MZ(object):
    """Numerical class to hold m/z values and do rich comparisons 
    and (some) math using the ppm error ranges."""
    
    __slots__ = ['m', 'z','mz', 'ppm', 'error', 'mz_range']
    
    def __init__(self,m,z=1,ppm=5):
        '''
        Can be instantiated using just m if z = 1 and ppm = 5.
        
        Args:
            m : mass 
            z : charge (defaults to 1 making mz = m)
            ppm : ppm error/precision of instrument (default = 5)

        '''
        self.m = m
        self.z = z
        self.ppm = ppm
        self.mz = m/z
        self.error = self.mz * (self.ppm / 10e6)
        self.mz_range = (
            (self.mz - self.error), 
            (self.mz + self.error))

    def __lt__(self, other):
        if isinstance(other,MZ):
            other = other.mz_range[0]

        if self.mz_range[1] < other:
            return True
        else:
            return False

    def __le__(self, other):
        if self.lt(other) or self.__eq__(other):
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

    def __str__(self):
        return "{mz:.4f} +/- {error:.4g}".format(mz=self.mz, error=self.error)

    def __repr__(self):
        return "{classn}(m={m},z={z},ppm={ppm})".format(
            classn = self.__class__.__name__,
            m = self.m,
            z = self.z,
            ppm = self.ppm)


