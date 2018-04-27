@cython.boundscheck(False)
@cython.wraparound(False)
def bin_search(np.ndarray[np.float_t, ndim=1] arr, float searchv):
    '''
    Binary search in cython.
    '''
    cdef int m
    cdef int minidx = 0
    cdef int maxidx = arr.shape[0] - 1
    while True:
        if maxidx < minidx:
            return arr.shape[0]
        m = (minidx + maxidx) // 2
        if arr[m] < searchv:
            minidx = m + 1
        elif arr[m] > searchv:
            maxidx = m - 1
        else:
            return m