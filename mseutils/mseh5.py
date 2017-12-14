'''
Work in progress...

This should be a mapping to HDF5 file format which I think is a natural fit for this MS data..
Might be overkill..
'''

import tables

from . import *

class MZ(tables.IsDescription):
    mz = tables.Float64Col()
    ppm = tables.Float64Col() 
    z = tables.Int8Col()
    intensity = tables.Float64Col()

